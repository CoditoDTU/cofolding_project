import pandas as pd
import numpy as np
import shutil
import zipfile
from pathlib import Path
from boltz_predictor import BOLTZPredictor

def process_batch_predictions(file_path: str):
    """
    Processes a CSV file for batch BOLTZ predictions.

    Args:
        file_path (str): The path to the input CSV file.
    """
    try:
        df = pd.read_csv(file_path)
        print(f"Loaded {len(df)} rows from {Path(file_path).name}")

        required_columns = ['Job_num', 'SMILES', 'Protein_sequence']
        if not all(col in df.columns for col in required_columns):
            raise KeyError(f"Missing required columns. Please ensure the CSV contains: {', '.join(required_columns)}")

        if 'IC50' not in df.columns:
            df['IC50'] = np.nan

        cif_temp_dir = Path("./temp_cifs")
        cif_temp_dir.mkdir(exist_ok=True)
        print(f"Temporary directory for CIF files created at {cif_temp_dir}")

        predictor = BOLTZPredictor()

        for index, row in df.iterrows():
            job_num = row['Job_num']
            protein_sequence_raw = row['Protein_sequence']
            ligand_smiles = row['SMILES']
            job_name = f"job_{job_num}"

            if isinstance(protein_sequence_raw, str):
                sequence_lines = protein_sequence_raw.split('\n')
                if len(sequence_lines) > 1 and sequence_lines[0].startswith('>'):
                    protein_sequence = "".join(sequence_lines[1:]).strip()
                else:
                    protein_sequence = protein_sequence_raw.strip()
            else:
                protein_sequence = ""

            if not protein_sequence or not ligand_smiles:
                print(f"Skipping row {index} (Job_num: {job_num}) due to missing data.")
                df.loc[index, 'IC50'] = "Skipped: Missing data"
                continue

            print(f"\nProcessing row {index} (Job_num: {job_num}):")
            print(f"Protein: {protein_sequence[:50]}... ({len(protein_sequence)} residues)")
            print(f"Ligand: {ligand_smiles}")

            input_data = {
                "version": 1,
                "sequences": [
                    {"protein": {"id": "A", "sequence": protein_sequence.upper().strip()}},
                    {"ligand": {"id": "B", "smiles": ligand_smiles.strip()}}
                ],
                "properties": [{"affinity": {"binder": "B"}}]
            }

            result = predictor.run_prediction(input_data, job_name)

            if result["success"]:
                metrics = predictor.extract_metrics(result["results"])
                if "ic50_nM" in metrics:
                    df.loc[index, 'IC50'] = metrics['ic50_nM']
                print(f"  Success! IC50 (nM): {metrics.get('ic50_nM', 'N/A')}")

                cif_files = result["files"].get("cif", [])
                for cif_file_path in cif_files:
                    original_cif_path = Path(cif_file_path)
                    new_cif_filename = f"{job_num}.cif"
                    target_cif_path = cif_temp_dir / new_cif_filename
                    try:
                        shutil.copy(original_cif_path, target_cif_path)
                        print(f"  Saved CIF file to: {target_cif_path}")
                        original_cif_path.unlink(missing_ok=True)
                    except Exception as e:
                        print(f"  Error managing CIF file {original_cif_path}: {e}")

                if "output_dir" in result and Path(result["output_dir"]).exists():
                    try:
                        shutil.rmtree(result["output_dir"])
                        print(f"  Cleaned up BOLTZ output directory: {result['output_dir']}")
                    except Exception as e:
                        print(f"  Error cleaning up directory {result['output_dir']}: {e}")
            else:
                error_message = result.get('error', 'Unknown error')
                print(f"  Prediction failed for row {index} (Job_num: {job_num}): {error_message}")
                df.loc[index, 'IC50'] = f"Error: {error_message}"

        file_name_stem = Path(file_path).stem
        output_csv_filename = f"{file_name_stem}_results.csv"
        output_csv_path = Path(f"./{output_csv_filename}")
        df.to_csv(output_csv_path, index=False)
        print(f"\nProcessing complete. Results CSV saved to: {output_csv_path}")

        final_zip_filename = f"{file_name_stem}_results_and_cifs.zip"
        final_zip_path = Path(f"./{final_zip_filename}")

        with zipfile.ZipFile(final_zip_path, 'w') as zipf:
            zipf.write(output_csv_path, arcname=output_csv_filename)
            for cif_file in cif_temp_dir.iterdir():
                if cif_file.suffix == '.cif':
                    zipf.write(cif_file, arcname=cif_file.name)

        print(f"Results CSV and CIF files zipped to: {final_zip_path}")

        shutil.rmtree(cif_temp_dir)
        output_csv_path.unlink(missing_ok=True)
        print("Cleaned up temporary files.")

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except KeyError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    # Replace 'your_file.csv' with the path to your CSV file
    input_file = 'your_file.csv'
    process_batch_predictions(input_file)