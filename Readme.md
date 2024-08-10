# Docking Pipeline
The Docking Pipeline offers an efficient and streamlined solution for molecular docking and MD simulations, leveraging Docker for easy setup. It includes step-by-step instructions for building and running the Docker container, executing Jupyter notebooks for data preprocessing, docking, and result extraction. The pipeline is optimized for high-performance computing, using parallel processing and GPU acceleration, and features advanced filtering steps to ensure only the most promising compounds are analyzed. 

![Docking PDF](devops/docking_image.png)

## Project Structure

- `Makefile`: Scripts to build and manage Docker containers.
- `devops`: Directory containing Docker configurations and additional resources.
  - `docker-compose.yml`: Docker Compose configuration file.
  - `jupyterlab`: Contains Dockerfile and environment configuration for JupyterLab.
  - `Filtration_flowchart_output.png`: Output image of the filtration flowchart.
  - `docking_image.png`: Image related to docking process.
  - `filtration_image.png`: Image related to filtration process.
- `main`: The main directory contains Jupyter notebooks for various processing steps and filtration operations, and additional subdirectories for specific tasks.
  - `1_process_data.ipynb`: Data processing notebook.
  - `2_unidock.ipynb`: UniDock processing notebook.
  - `3_check_progress.ipynb`: Progress checking notebook.
  - `4_extract_data.ipynb`: Data extraction notebook.
  - `5_frag_filteration.ipynb`: Fragment filtration notebook.
  - `6_smol_filteration.ipynb`: Small molecule filtration notebook.
  - `7_protac_filteration.ipynb`: PROTAC filtration notebook.
  - `frag_filteration`: Directory for fragment filtration resources.
  - `protac_filteration`: Directory for PROTAC filtration resources.
  - `smol_filteration`: Directory for small molecule filtration resources.
  - `test`: Directory for test scripts and resources.
  - `z_docking_files`: Additional docking scripts and resources.
  - `z_filteration_files`: Additional filtration scripts and resources.


## Docker Setup and Usage
Note: Make sure you have Docker installed on your system.
1. ![#ffa500](https://via.placeholder.com/15/ffa500/ffa500.png) Clone the repository to your local machine: `git clone git@github.com:KailashLohar/Docking_and_Filteration_Pipeline.git`

2. ![#ffa500](https://via.placeholder.com/15/ffa500/ffa500.png) Navigate to the project directory "**Docking_and_Filteration_Pipeline**" and run the following commands to build the Docker image :`make build-image`

3. ![#ffa500](https://via.placeholder.com/15/ffa500/ffa500.png) Start the container using command: `make start-container`

4. ![#ffa500](https://via.placeholder.com/15/ffa500/ffa500.png) Enter the Container using command: `make enter-container`

5. ![#ffa500](https://via.placeholder.com/15/ffa500/ffa500.png) Stop the Container using command: `make stop-container`



## Pipeline Instruction
Once you are inside Jupyterlab container, you have to execute four notebooks to get the final outputs suitable for Plip and MD simulations, "1_process_data", "2_unidock", "3_check_progress" and "4_extract_data" simultaneously.

### Instruction for running "1_process_data":

1. ![#DA70D6](https://via.placeholder.com/15/DA70D6/DA70D6.png) Create a new folder and name it for e.g `my_folder`.
   
2. ![#DA70D6](https://via.placeholder.com/15/DA70D6/DA70D6.png) Place your CSV file (e.g., **my_csv.csv**) inside the `my_folder`. Make sure it has atleast 5 compounds
   
3. ![#DA70D6](https://via.placeholder.com/15/DA70D6/DA70D6.png) In cell No. 2 of the notebook, mention your `my_folder` in **sub_folder_name** and csv file name in **file_name**.
4. ![#DA70D6](https://via.placeholder.com/15/DA70D6/DA70D6.png) In cell No. 3, examine your DataFrame to identify the column containing compound names and the second column containing SMILES.
7. ![#DA70D6](https://via.placeholder.com/15/DA70D6/DA70D6.png) In cell No. 4, rename these columns to **Name** and **SMILES** respectively.
8. ![#DA70D6](https://via.placeholder.com/15/DA70D6/DA70D6.png) After running the code, you will find a new file named **my_csv_out.csv** in the `my_folder`, which will be suitable for Docking.


### Step by step instruction for running "2_unidock":
1. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `Load Data`: Place following files inside `my_folder`
    - ligands csv file (e.g., ligands_smiles_out.csv, mention the file that is obtained from previous preprocess notebook)
    - protein PDBQT file (e.g., **my_pdbqt.pdbqt**)
    - protein PDB file (e.g., **my_pdb.pdb**)
    - configuration file (e.g., **my_conf.txt**) inside the `my_folder`.
    - threshold (Threshold can be static or dynamic.)
    - factor (factor can be 1, 2 or 3.) (Static means, if you pass value let say -7 (numeric) then all ligands having BA less than -7 will be going for posebusters filteration and here any factor value become unusefull or of no value. Dynamic means, if you set threshold = 'dynamic' and factor = 2 then distribution of BA will be calculated and all ligands having BA less than μ - 2σ will be pass for posebusters filteration.)

2. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `SMILES to SDF`: The code reads the SMILES strings from the DataFrame, generates molecular structures with a given number of conformations, and saves these structures as SDF files in a designated output directory (pipeline_files/1_sdf). The process is parallelized using multiple CPU cores to speed up the conversion. Each row in the DataFrame is processed to convert the SMILES to a 3D structure, and the resulting SDF files are named based on the Name column from the DataFrame.

3. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `SDF to Mol2`: This step runs a Bash script (1_sdf_to_mol2.sh) which performs the following steps: it minimizes the structures in batches of 20,000 using the obminimize tool, then converts the minimized SDF files to MOL2 format in batches of 20,000 using obabel. The output MOL2 files are saved in a designated folder (pipeline_files/2_mol2). The process ensures efficient handling of large numbers of files by processing them in manageable batches.

4. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `Format Mol2 files`:  This step reads each MOL2 file, adds sequential atom numbers to each atom type for clarity, and saves the formatted files to a new directory. The process involves identifying atom sections in the MOL2 files, assigning unique numbers to each atom type, and maintaining proper alignment. The formatted MOL2 files are saved in a designated output folder (pipeline_files/2_mol2_format) using parallel processing to handle multiple files efficiently.

5. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `Mol2 to PDBQT`: This step runs a Bash script (2_mol2_to_pdbqt.sh) which performs the following steps: it converts the formatted MOL2 files to PDBQT format in batches of 20,000 using obabel. The output PDBQT files are saved in a designated folder (pipeline_files/3_pdbqt). The process ensures efficient handling of large numbers of files by processing them in manageable batches, and the completion of the conversion is confirmed by a printed message.

6. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `PDBQT Filter`: There are threes substeps in PDBQT Filter steps.
    - `PDBQT to SMILES`: This step runs a Bash script (3_pdbqt_to_smiles.sh) which performs the following steps: it converts the PDBQT files to SMILES format in batches of 20,000 using obabel. The output SMILES files are saved in a designated folder (pipeline_files/4_smiles). The process ensures efficient handling of large numbers of files by processing them in manageable batches, and the completion of the conversion is confirmed by a printed message.
    - `Filter Compounds`: The code reads SMILES files and merges them with an input CSV to create a DataFrame, calculates Dice similarity scores between SMILES strings, filters out compounds with a perfect similarity score, and saves the filtered results to pipeline_files/1_compounds_for_docking.csv.
    - `Copy Correct PDBQT Files`: The code copies the corresponding PDBQT files of the filtered compounds to a pipeline_files/5_pdbqt_for_docking folder, ensuring only the correct PDBQT files are prepared for docking.

7. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `UniDock Docking`: There are two substeps involved in this step.
    - `Batch Creation`: The code creates batches of PDBQT files from a specified folder (pipeline_files/5_pdbqt_for_docking). It reads all PDBQT files, splits them into batches of 20,000 files each, and writes the paths of each batch to a separate text file (unidock_pdbqt_batch_X.txt).
    - `Docking Execution`: The code iterates over each batch file, constructs and runs a UniDock command to perform docking using the specified receptor and configuration files. The docking results are saved in a designated output directory (pipeline_files/6_pdbqt_out), and logs are maintained for each batch (unidock_output_batch_X.txt). The process ensures all batches are docked and results are saved accordingly.

8. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `Extract Affinity values`: This step reads PDBQT files from the pipeline_files/6_pdbqt_out directory, extracts affinity values from each file, and compiles these values along with the corresponding compound names into a list. It then writes this list to a CSV file (pipeline_files/2_extract_affinity_from_pdbqt.csv). The process ensures that all affinity values are extracted and saved in a consolidated CSV file for further analysis.

9. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `Extract Compounds based on Affinity threshold`: There are two steps here.
    - `Filter Compounds Based on Affinity`: The code reads affinity scores from a CSV file, calculates a dynamic or static threshold, filters compounds with affinity values below this threshold, and copies the corresponding PDBQT files to a new directory (pipeline_files/7_pdbqt_out_threshold).
    - `Extract Model 1 Content`: The code reads the filtered PDBQT files, extracts content up to "ENDMDL", and saves this extracted content to another directory (pipeline_files/8_pdbqt_out_threshold_m1), ensuring only the relevant parts of the PDBQT files are retained for further processing.

10. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `PDBQT to Mol2 and Mol2 to SDF`: This step converts PDBQT files from the pipeline_files/8_pdbqt_out_threshold_m1 directory to MOL2 format and then to SDF format, processing files in batches of 20,000. It creates two output directories (pipeline_files/9_mol2_out and pipeline_files/9_sdf_out), and uses the obabel tool to perform the conversions in batches, ensuring efficient handling of large numbers of files. The conversion process is confirmed by a printed message indicating the completion and the location of the saved SDF files.

11. ![#c5f015](https://via.placeholder.com/15/c5f015/c5f015.png) `PoseBusters Filter`: Two steps are there.
    - `PoseBuster filter`: The code runs a Bash script that filters SDF files using PoseBusters, checking for docking compatibility with a specified PDB file. It reads SDF files from the pipeline_files/9_sdf_out directory, processes them in parallel, and saves the filtration results in a CSV file (pipeline_files/4_pb_out.csv). If the output CSV file already exists, it is removed before starting the process. The progress and completion of the PoseBusters filtration are confirmed by printed messages.
    - `Process PoseBusters Results`: The step reads the PoseBusters results from a CSV file (4_pb_out.csv), drops specific columns and rows, renames columns, calculates the number of passes for each compound, and saves the cleaned results to a new CSV file (5_pb_out.csv).
    - `Generate Final Output`: The code merges the processed PoseBusters data with affinity scores and input SMILES data, filters compounds based on the number of passes, calculates the efficiency of each compound, and saves the final results to an output CSV file (output.csv). This ensures that only the most efficient compounds passing the PoseBusters criteria are included in the final output.



### Instruction for running "3_check_progress":

1. ![#4169E1](https://via.placeholder.com/15/367588/367588.png) This notebook monitors the progress of various file conversions and filtrations in a pipeline. It checks the number of files at each stage: from SMILES to SDF, SDF to MOL2, MOL2 to PDBQT, and PDBQT to SMILES, including the count of PDBQT files filtered for docking. It also tracks the number of files post-docking, after applying a threshold filter, and conversions from PDBQT to MOL2 and MOL2 to SDF. Additionally, it monitors the progress of PoseBusters filtration and the total number of final output compounds, ensuring each step's completion and correctness.



### Explanation of "4_extract_data" notebook:

1. ![#FFDB58](https://via.placeholder.com/15/FFDB58/FFDB58.png) `Extract ligands in pdbqt format for Analysis`: Mention the folder name and csv path from where you want to extract data and also the number of files to extract. This section copies and zips a specified number of PDBQT files based on the top entries in a given CSV file. It reads the CSV to get the names of the top compounds, copies their corresponding PDBQT files from the pipeline_files/7_pdbqt_out_threshold directory to a new directory (pdbqt_extracted), and then zips these files into a single archive. After creating the ZIP file, the temporary directory used for copying the PDBQT files is deleted. The process ensures that the top num_files PDBQT files are efficiently packaged into a ZIP archive for easy distribution or further analysis.

2. ![#FFDB58](https://via.placeholder.com/15/FFDB58/FFDB58.png) `Extract protein + ligand in pdb format for PLIP Analysis`: Mention the folder name and csv path from where you want to extract data and also the number of files to extract. The code converts PDBQT files to PDB format and then zips a specified number of these converted files along with a protein structure file. It reads a CSV file to get the top compounds, converts their corresponding PDBQT files from pipeline_files/8_pdbqt_out_threshold_m1 to PDB format, and appends the protein's ATOM lines to these PDB files. The converted files are temporarily stored in a new directory (pdbqt_extracted). Finally, the code creates a ZIP archive containing these PDB files and deletes the temporary directory. This process ensures that the top num_files PDB files, along with the protein structure, are efficiently packaged into a ZIP archive for easy distribution or further analysis.

### Execution time
Virtual screening will take approximately 26 hrs for 1 million compounds if a system is having 32 cores CPU and Nvidia GeForce RTX 4090 GPU.

### Built With
- NVIDIA CUDA_12.2 - Parallel computing platform and programming model
- Miniconda - Python package management
- Boost - C++ libraries
- NVIDIA Container Runtime - NVIDIA container runtime

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Filteration Pipeline
The Filtration Pipeline is a comprehensive solution designed to streamline the filtering of molecular compounds through distinct processes tailored for FRAG, SMOL, and PROTAC. This pipeline features steps such as loading and cleaning data, removing salts and duplicates, applying molecular descriptor and catalog filters, and eliminating unwanted ring structures. Each step meticulously refines the dataset, ensuring only the most relevant compounds are retained. Visual flowcharts enhance understanding by mapping out the entire filtration workflow. This pipeline is essential for efficient and effective molecular compound analysis, providing clear, step-by-step instructions to achieve high-quality results.

![Filteration PNG](devops/filtration_image.png)

## Pipeline Instruction
In filtration, there are three separate pipelines for FRAG, SMOL, and PROTAC. All these pipelines are identical except for step 4: Molecular Descriptor Filtration.

### Step by step explanation of "5_frag_filteration" notebook:
1. ![#00FFFF](https://via.placeholder.com/15/00FFFF/00FFFF.png) `Load Data`: The code reads a CSV file from the frag_filteration folder, which contains compound names and SMILES strings. It then selects only the columns with the compound names and SMILES strings, renames these columns to 'Name' and 'SMILES' respectively, and displays the first few rows and the shape of the DataFrame before and after the renaming. This process ensures that only the relevant columns are retained and properly labeled for further analysis.

2. ![#00FFFF](https://via.placeholder.com/15/00FFFF/00FFFF.png) `Remove Salts`: This step processes a DataFrame to remove salts from SMILES strings. It converts each SMILES string to a molecule, breaks bonds, neutralizes it, and filters out nonorganic and salt fragments. The cleaned fragments are then recombined into a new SMILES string, standardized, and any errors during this process are handled.

3. ![#00FFFF](https://via.placeholder.com/15/00FFFF/00FFFF.png) `Remove Duplicates`: The code defines a function process_smiles_dataframe to clean a DataFrame by canonicalizing SMILES strings and removing duplicates based on canonical SMILES and InChI keys. It deep copies the DataFrame, resets its index, canonicalizes SMILES, and removes duplicates. The function prints the number of removed and remaining data points, and returns the cleaned DataFrame and the final count. This function is then applied to df_remove_salts, storing the cleaned DataFrame and counts in df_duplicates and after_duplicates_count.

4. ![#00FFFF](https://via.placeholder.com/15/00FFFF/00FFFF.png) `Molecular Descriptor Filteration`: The code defines a function frag_calculate_and_filter_properties to calculate molecular properties for each SMILES string and filter compounds based on defined thresholds. It computes properties like molecular weight, cLogP, HBA, HBD, and NRotB, then filters compounds meeting specific criteria. The function prints the number of removed and remaining data points, returning the filtered DataFrame and the final count. This function is applied to df_duplicates, resulting in df_filtered and after_filter_properties.

5. ![#00FFFF](https://via.placeholder.com/15/00FFFF/00FFFF.png) `RDKit Catalog Filteration`: This step defines a function catalog_filter to filter compounds based on various catalog filters such as PAINS, NIH, ZINC, and others. It checks each SMILES string against these catalogs, removes compounds with matches, and returns the filtered DataFrame and the final count. The function prints the number of removed and remaining data points, and saves the intermediate results to a CSV file. This function is applied to df_filtered, resulting in df_catlog and after_catlog_filter.

6. ![#00FFFF](https://via.placeholder.com/15/00FFFF/00FFFF.png) `Remove rings with consecutive heteroatoms`: If the flag execute_remove_rings_with_consecutive_heteroatoms is set to true, the function remove_rings_with_consecutive_heteroatoms is applied to df_catlog to filter out compounds containing rings with three consecutive heteroatoms. The function checks each SMILES string, removes compounds with such rings, and returns the filtered DataFrame and the final count. If the flag is false, the original DataFrame df_catlog is used without additional filtering. The function prints the number of removed and remaining data points, resulting in df_con_het and after_consecutive_het_atoms.

7. ![#00FFFF](https://via.placeholder.com/15/00FFFF/00FFFF.png) `Remove max three conjugated fused rings`: The code defines a function remove_max_three_conjugated_fused_rings that filters out compounds with more than three conjugated fused rings. It checks each SMILES string, identifies fused rings with double bonds, and retains compounds with three or fewer such rings. The function prints the number of removed and remaining data points, returning the filtered DataFrame and the final count. This function is applied to df_con_het, resulting in df_conj_rings and after_conjugated_fused.

8. ![#00FFFF](https://via.placeholder.com/15/00FFFF/00FFFF.png) `Remove rings with more than seven members`: The code defines a function remove_rings_with_more_than_seven_members that filters out compounds containing rings with more than seven members. It checks each SMILES string for such rings and removes the corresponding compounds. The filtered DataFrame is saved to output_filteration.csv, and the function prints the number of removed and remaining data points. This function is applied to df_conj_rings, resulting in df_final and after_remove_seven_member.

9. ![#00FFFF](https://via.placeholder.com/15/00FFFF/00FFFF.png) `Filteration Flowchart`: The code defines a function create_flowchart to visualize the filtration process. It calculates the reduction in compound counts at each step, creates a flowchart using graphviz, and labels each node with the step name and compound count. The flowchart includes steps like removing salts, duplicates, filtering by molecular descriptors, catalog filtration, and removing specific ring structures. The final flowchart is saved as a PNG file in the specified folder, summarizing the filtration workflow.

![Filtration Flowchart](devops/Filtration_flowchart_output.png)


## Author

- Kailash Lohar  <kailashlohariitkgp@gmail.com>


## Reference
Yu, Y., Cai, C., Wang, J., Bo, Z., Zhu, Z., & Zheng, H. (2023). Uni-Dock: GPU-Accelerated Docking Enables Ultralarge Virtual Screening. Journal of Chemical Theory and Computation. https://doi.org/10.1021/acs.jctc.2c01145

## License

This project is licensed under the MIT License - see the LICENSE file for details.
