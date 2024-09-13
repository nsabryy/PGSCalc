## Config File Description

This configuration file is used to control the behavior of the PGSCalc application. Below is a description of each section and field in the configuration file:

### `VCFQuery`
This section contains settings related to the querying of VCF (Variant Call Format) files.

- **`input_dir`**: 
  - *Description*: The directory where the input VCF files are located.
  - *Type*: String containing filepath.
  - *Example*: `/N/project/biobank/imputed-from-dongbing/reheader/`

- **`remove_imputed`**: 
  - *Description*: A flag to indicate whether imputed variants should be removed.
  - *Type*: Boolean

- **`remove_filtered`**: 
  - *Description*: A flag to indicate whether filtered variants should be removed.
  - *Type*: Boolean

- **`sample_ids`**: 
  - *Description*: A list of sample IDs to be queried.
  - *Type*: List of strings OR .txt file
  - *Example*: `["23514","41689"]`, `path/to/samples.txt`

- **`override`**: 
  - *Description*: A flag to indicate whether previous results should be overridden. Set this to `true` if there have been any changes to your input files or PGS scoring files. Additionally, if computation is interrupted, these objects may be corrupted, so you want to set the override to `true` in that case as well.
  - *Type*: Boolean

### `PGScore`
This section contains settings related to the polygenic score (PGS) files.

- **`pgs_ids`**: 
  - *Description*: A list of PGS IDs to be used for scoring.
  - *Type*: List of strings
  - *Example*: `["PGS000001", "PGS000008", "PGS000035"]`

- **`genome_build`**: 
  - *Description*: The genome build version to be used for scoring download, this must match your input data. 
  - *Type*: String
  - *Example*: `GRCh37`

- **`override`**: 
  - *Description*: A flag to indicate whether previous PGS scoring results should be overridden. Similar to the `VCFQuery` override, set this to `true` if there have been any changes to the PGS scoring files or if there has been an interruption in a past computation.
  - *Type*: Boolean
  - *Example*: `false`

**Note**: You may use your own scoring file for computation rather than a score file from the PGSCatalog. To do this, add `scoring_file_directory` to `PGScore`. You may delete the above fields; even if they are filled out, if this field is included *and* has a value, the pipeline will automatically use the local scorefile and ignore the other fields.

### `PGSCalc`
This section contains settings related to the calculation and output of the polygenic scores.

- **`parallel`**: 
  - *Description*: A flag to indicate whether multitprocessing should be used.
  - *Type*: Boolean
  - *Example*: `true`

- **`output_dir`**: 
  - *Description*: The directory where the output results will be saved.
  - *Type*: String
  - *Example*: `./scoring_results/testing/`

- **`summary_statistics`**: 
  - *Description*: A flag to indicate whether summary statistics should be calculated and saved.
  - *Type*: Boolean
  - *Example*: `true`

- **`include_snps_per_sample`**:
  - *Description*: Determines whether you would like to include the SNPs (stored as (chrom, pos, score)) found for each participant in the raw scoring file. Helpful to compare cohorts. Note that for large scoring files, this will massively increase the size of your raw scoring output so should only be used in a testing environment.
  - *Type*: Boolean



### Example Configuration File
```json
{
    "VCFQuery": {
        "input_dir": "/N/project/biobank/imputed-from-dongbing/reheader/",
        "remove_imputed": false,
        "remove_filtered": false,
        "sample_ids": ["23514","41689"],
        "override": false
    },
    "PGScore": {
        "pgs_ids": ["PGS000001", "PGS000008", "PGS000035"],
        "genome_build": "GRCh37",
        "override": false
    },
    "PGSCalc": {
        "multithreading": true,
        "output_dir": "./scoring_results/testing/",
        "input_dir": "",
        "summary_statistics": true
    }
}
```
### Summary Statistics

- **Summary Report**: After processing each PGS score, a summary report is generated.
  - **Match Percentage**: Percentage of alleles found compared to total alleles queried.
  - **Top 10% Scores**: Threshold for the top 10% scores and the list of samples in this category.
  - **Score Distribution**: A histogram plot showing the distribution of scores.

The summary report includes:
- A text file with match percentage, top 10% score threshold, and top 10% samples.
- A histogram plot showing the distribution of scores.

## Output Filenames and Directory Structure

The output filenames and directory structure are organized to ensure clarity and consistency. Below is the description of the output format:

### Output Filename Format

The output filenames follow the structure:

`{pgs_id}_{normalized_input_dir}_{optional_uid_for_subsample}_{1 if imputed removed, 0 otherwise}_{1 if filtered removed, 0 otherwise}`


- **`pgs_id`**: The polygenic score ID.
- **`normalized_input_dir`**: The normalized input directory name, where special characters are replaced with underscores.
- **`optional_uid_for_subsample`**: An optional unique identifier for subsamples, if applicable.
- **`1 if imputed removed, 0 otherwise`**: Indicates whether imputed variants were removed (1) or not (0).
- **`1 if filtered removed, 0 otherwise`**: Indicates whether filtered variants were removed (1) or not (0).

### Directory Structure

The directory structure for the outputs is as follows:

1. **Top-level directory**: The normalized input directory name within the `PGSCalc/scoring_results` directory.
2. **Subdirectory**: The polygenic score ID (`pgs_id`).
3. **Filename**: The formatted name as described above.

### Example Directory Structure
```
PGSCalc/scoring_results
└─── input_directory-normalized_0_0
    └─── PGS000001
        │   PGS000001_N_project_biobank_imputed-from-dongbing_reheader_0_0.tsv
        │   PGS000001_N_project_biobank_imputed-from-dongbing_reheader_0_0_summary.txt
        │   PGS000001_N_project_biobank_imputed-from-dongbing_reheader_0_0_score_distribution.png
```

The output location is also included at the end of the log file. 

**Warning** if you run scoring multiple times on the same inputs and scoring file, this will overwrite previous results unless you change the `output_directory` in `PGSCalc` in the config file.

## TODO
- scoring files from multiple identifiers (trait, publication)
- checking distribution ()