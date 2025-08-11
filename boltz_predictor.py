import subprocess
import json
import tempfile
import yaml
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np


class BOLTZPredictor:
    """Main interface for BOLTZ 2 predictions"""

    def __init__(self, output_dir: str = "output"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

    def run_prediction(self, input_data: dict, job_name: str,
                      use_msa: bool = True, verbose: bool = False) -> dict:
        """
        Execute BOLTZ prediction

        Args:
            input_data: YAML-compatible input dictionary
            job_name: Unique identifier for this job
            use_msa: Whether to use MSA server
            verbose: Print detailed output

        Returns:
            Dictionary containing results and file paths
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(input_data, f)
            yaml_path = f.name

        output_path = self.output_dir / job_name
        cmd = ["boltz", "predict", yaml_path, "--out_dir", str(output_path)]

        if use_msa:
            cmd.append("--use_msa_server")

        try:
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                return {"success": False, "error": stderr or stdout}

            result_dirs = list(output_path.glob("boltz_results_*"))

            files = {"pdb": [], "cif": [], "json": []}
            for result_dir in result_dirs:
                files["pdb"].extend(list(result_dir.rglob("*.pdb")))
                files["cif"].extend(list(result_dir.rglob("*.cif")))
                files["json"].extend(list(result_dir.rglob("*.json")))

            results = {}
            for json_file in files["json"]:
                try:
                    with open(json_file) as f:
                        results[json_file.name] = json.load(f)
                except Exception as e:
                    if verbose:
                        print(f"Warning: Could not read {json_file.name}: {e}")

            return {
                "success": True,
                "output_dir": str(output_path),
                "files": files,
                "results": results,
                "stdout": stdout if verbose else None
            }

        except Exception as e:
            return {"success": False, "error": str(e)}
        finally:
            Path(yaml_path).unlink(missing_ok=True)

    def extract_metrics(self, results: dict) -> dict:
        """Extract key metrics from BOLTZ results"""
        metrics = {}

        for filename, data in results.items():
            if "confidence" in filename and isinstance(data, dict):
                metrics.update({
                    "confidence_score": data.get("confidence_score", 0),
                    "ptm": data.get("ptm", 0),
                    "iptm": data.get("iptm", 0),
                    "plddt": data.get("complex_plddt", 0)
                })
            elif "affinity" in filename and isinstance(data, dict):
                log_ic50_uM = data.get("affinity_pred_value", 0)
                ic50_uM = 10 ** log_ic50_uM
                ic50_nM = ic50_uM * 1000
                ic50_M = ic50_uM * 1e-6
                pic50 = -np.log10(ic50_M) if ic50_M > 0 else 0
                delta_g_kcal = (6 - log_ic50_uM) * 1.364
                kd_uM = ic50_uM / 2
                kd_nM = kd_uM * 1000
                kd_M = kd_uM * 1e-6
                pkd = -np.log10(kd_M) if kd_M > 0 else 0

                metrics.update({
                    "boltz_affinity_value": log_ic50_uM,
                    "log_ic50_uM": log_ic50_uM,
                    "ic50_uM": ic50_uM,
                    "ic50_nM": ic50_nM,
                    "pic50": pic50,
                    "kd_uM": kd_uM,
                    "kd_nM": kd_nM,
                    "pkd": pkd,
                    "delta_g_kcal": delta_g_kcal,
                    "affinity_prob": data.get("affinity_probability_binary", 0)
                })

        return metrics