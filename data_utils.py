import pandas as pd
import itertools
import os
from typing import Dict, List, Any



def create_combination_dataset(
    param_combinations: Dict[str, List[Any]],
    save: bool = False,
    output_dir: str = None,
    filename: str = None,
    add_id_column: bool = True
) -> pd.DataFrame:
    """
    Generates a Pandas DataFrame containing all unique combinations of parameters.

    This function takes a dictionary where keys are parameter names and values are
    lists of options for those parameters. It then computes the Cartesian product
    to generate a DataFrame of all possible unique combinations.

    Args:
        param_combinations (Dict[str, List[Any]]): A dictionary where keys
            are the parameter names (column headers) and values are lists
            of the parameter options.
        save (bool, optional): If True, the DataFrame will be saved to a CSV file.
            Defaults to False.
        output_dir (str, optional): The directory path where the file will be saved.
            Required if save is True.
        filename (str, optional): The name for the output CSV file.
            Required if save is True.
        add_id_column (bool, optional): If True, an 'id' column will be added
            to the DataFrame. Defaults to True.

    Returns:
        pd.DataFrame: A DataFrame where each row is a unique combination of
            the input parameters.

    Raises:
        ValueError: If `save` is True but `output_dir` or `filename` is not provided.
    """
    if save and (not output_dir or not filename):
        raise ValueError("If 'save' is True, both 'output_dir' and 'filename' must be provided.")

    # Get parameter names and their corresponding value lists
    param_names = list(param_combinations.keys())
    param_values = list(param_combinations.values())
    print(param_names)
    print(param_values)
    # Generate all unique combinations using itertools.product
    combinations = list(itertools.product(*param_values))
    total_combinations = len(combinations)
    print(f"Generated {total_combinations} unique combinations.")

    # Create the DataFrame
    df = pd.DataFrame(combinations, columns=param_names)

    # Add a unique ID column at the beginning if requested
    if add_id_column:
        df.insert(0, 'Job_num', range(1, total_combinations + 1))

    # Save the DataFrame to a CSV file if requested


    if save:
        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        csv_path = os.path.join(output_dir, filename)
        df.to_csv(csv_path, index=False)
        print(f"Successfully saved DataFrame to {csv_path}")

    return df



def merge_and_prepare_data(combinatorial_df: pd.DataFrame, 
                           homologs_df: pd.DataFrame, 
                           substrates_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merges combinatorial, homolog, and substrate dataframes and formats the result.

    This function performs two consecutive left merges to combine the three input 
    dataframes. It then renames the relevant columns and selects a final subset 
    of columns required for further analysis.

    Args:
        combinatorial_df (pd.DataFrame): The primary dataframe containing combinations of 
                                          substrates and homologs. It must include the columns 
                                          'substrates', 'homologs', and 'Job_num'.
        homologs_df (pd.DataFrame): The dataframe containing details about the homologs. 
                                    It must include the 'homolog' and 'aasequence' columns.
        substrates_df (pd.DataFrame): The dataframe containing details about the substrates. 
                                      It must include the 'substrate' and 'SMILES_charged' columns.

    Returns:
        pd.DataFrame: A new, clean dataframe with the merged and formatted data. The output 
                      contains the columns 'Job_num', 'SMILES', and 'Protein_sequence'.
    """
    # Merge the combinatorial dataframe with the substrates dataframe
    merged_df = pd.merge(
        combinatorial_df,
        substrates_df,
        left_on='substrates',
        right_on='substrate',
        how='left'
    )

    # Merge the result with the homologs dataframe
    final_df = pd.merge(
        merged_df,
        homologs_df,
        left_on='homologs',
        right_on='homolog',
        how='left'
    )

    # Rename the columns for clarity and consistency
    final_df = final_df.rename(columns={
        'SMILES_charged': 'SMILES',
        'aasequence': 'Protein_sequence'
    })

    # Select and reorder the final columns needed for the output
    final_df = final_df[[
        'Job_num',
        'homolog_id',
        'substrate_id'
        'SMILES',
        'Protein_sequence'
    ]]

    return final_df