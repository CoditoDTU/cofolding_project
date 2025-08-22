# Enzyme-Ligand Cofolding with Boltz-2

This repository provides a user-friendly Jupyter Notebook to perform cofolding of enzyme structures with their respective ligands and cofactors using the Boltz-2 model.

## Overview

The `boltz_pipeline.ipynb` notebook in this project guides you through the entire process of setting up and running a cofolding simulation. It simplifies the workflow by bundling all the necessary steps, from environment setup to executing the folding model. This tool is designed for researchers who want to model complex protein structures without getting bogged down in complex command-line operations.

---

## Prerequisites

Before you begin, ensure you have the following installed on your system:

* **Conda**: This project relies on a Conda environment to manage dependencies. We recommend installing **Miniforge** or **Anaconda**. You can find installation instructions for Miniforge [here](https://github.com/conda-forge/miniforge).
* **NVIDIA GPU**: Running Boltz-2 is computationally intensive. An NVIDIA GPU with CUDA installed is **highly recommended** for reasonable performance. The process may be extremely slow or fail on a CPU.

---

## Step-by-Step Guide

Follow these steps to get the notebook up and running on your local machine.

### 1. Clone the Repository

First, clone this repository to your local machine and navigate into the project directory.

```
git clone https://github.com/CoditoDTU/cofolding_project
cd cofolding_project
```

### 2. Set Up the Conda Environment
Next, create the Conda environment using the provided ```requirements.yaml``` file. This file contains all the necessary dependencies. Activate the environment once it's created.

```bash
# Create the environment from the file
conda env create -f requirements.yaml
# Activate the newly created environment
conda activate boltz2_env
```
### 3. Run the Notebook

With the `boltz2_env` environment active, you can now launch Jupyter and open the main pipeline notebook.

Once Jupyter opens in your browser, click on the file named boltz_pipeline.ipynb and follow the instructions inside.

Happy Folding! 🙏