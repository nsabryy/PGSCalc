# mean_rscore = sum(n=0..k){ abs(pop_mean_n - total_mean) }
# # stdv_rscore = sum(n=0..k){ abs(pop_stdv_n - total_stdv) }
# For the PGS, can you calculate the following:
# 1) The mean value for all participants (all populations)
# 2) For each population, the mean value.

# The final “robustness/transferability” value will be the sum of the differences between each population and the total mean. (Absolute value the difference).

# robustness = sum(n=0..k){ abs(pop_mean_n - total_mean) }
import numpy as np
import subprocess
import os

def find_tsv_files(directory):
    tsv_files = []
    for root, dirs, files in os.walk(directory):
        tsv_in_subdir = [os.path.join(root, file) for file in files if file.endswith('.tsv')]
        if len(tsv_in_subdir) == 1:
            tsv_files.append(tsv_in_subdir[0])
        elif len(tsv_in_subdir) > 1:
            print(f"Warning: Multiple .tsv files found in {root}. Skipping this subdirectory.")
        elif len(tsv_in_subdir) == 0:
            print(f"Warning: No .tsv file found in {root}. Skipping this subdirectory.")
    return tsv_files

def load_population_info(index_fp, superpopulation=True):
    superpopulations = {
        "CHB": "EAS", "JPT": "EAS", "CHS": "EAS", "CDX": "EAS", "KHV": "EAS",
        "CEU": "EUR", "TSI": "EUR", "FIN": "EUR", "GBR": "EUR", "IBS": "EUR",
        "YRI": "AFR", "LWK": "AFR", "MAG": "AFR", "MSL": "AFR", "ESN": "AFR",
        "ASW": "AFR", "ACB": "AFR", "GWD": "AFR",
        "MXL": "AMR", "PUR": "AMR", "CLM": "AMR", "PEL": "AMR",
        "GIH": "SAS", "PJL": "SAS", "BEB": "SAS", "STU": "SAS", "ITU": "SAS"
    }

    # Load population info for samples as needed
    sample_to_population = {}
    with open(index_fp, 'r') as f:
        # Skip metadata lines starting with ##
        line = f.readline().strip()
        while line.startswith("##"):
            line = f.readline().strip()
        
        # Parse the header
        if line.startswith("#"):
            header = line.lstrip('#').strip().split('\t')
            pop_idx = header.index("POPULATION")
            sample_id_idx = header.index("SAMPLE_NAME")
        else:
            raise ValueError("Index file header not found.")
        
        # Load relevant population data
        for line in f:
            fields = line.strip().split('\t')
            sample_id = fields[sample_id_idx]
            population = fields[pop_idx]
            if superpopulation:

                sample_to_population[sample_id] = superpopulations[population]
            else:
                sample_to_population[sample_id] = population
    
    return sample_to_population

def get_population_distributions(scoring_results_fp, index_fp):
    # Load the sample-to-population mapping
    sample_to_population = load_population_info(index_fp)
    
    # Dictionary to store population-specific score lists
    pop_distributions = {'ALL':[]}
    all_scores = []
    with open(scoring_results_fp, 'r') as f:
        # Skip the header line in the scoring file
        header = f.readline()
        for line in f:
            fields = line.strip().split('\t')
            sample_id = fields[1]
            score = float(fields[3])
            
            population = sample_to_population.get(sample_id, "UNKNOWN")

            if population not in pop_distributions:
                pop_distributions[population] = []
            pop_distributions[population].append(score)
            pop_distributions['ALL'].append(score)

    return pop_distributions

def calculate_mean_robustness(pop_distributions):
    total_mean = np.average(pop_distributions.pop('ALL'))
    robustness = 0
    for _, scores in pop_distributions.items():
        robustness += abs(np.average(scores)-total_mean)
    return robustness

def calculate_stdv_robustness(pop_distributions):
    # TODO
    pass 

def main():
    base_dir = '/N/project/compgen/PGSCalc/scoring_results/9e09114e_0_0_union/EFO_0004736'
    input_files = find_tsv_files(base_dir)
    index_fp = '/N/project/compgen/shared/resources/1000_genomes/20220422_3202_phased_SNV_INDEL_SV/1000genomes.sequence.index'
    total_files = len(input_files)
    score_robustness = {}

    for i, input_fp in enumerate(input_files, start=1):
        print(f"[{i}/{total_files}] Processing {input_fp}...", end='', flush=True)
        pop_distributions = get_population_distributions(input_fp, index_fp)
        robustness = calculate_mean_robustness(pop_distributions)
        if input_fp in score_robustness:
            print("ERROR! FP ALREADY PROCESSED. NOT SAVING NEW SCORE.")
        score_robustness[input_fp] = robustness
        print(f"\r[{i}/{total_files}] Processing {input_fp}... Done.")
    
    output_base = 'robustness.txt'
    output_fp = os.path.join(base_dir, output_base)

    with open(output_fp, 'w') as f:
        f.write(f"Base Directory: {base_dir}\n\n")
        f.write("filepath\tmean_robustness\n")

        for key, value in score_robustness.items():
            stripped_path = key.replace(base_dir + '/', '')
            f.write(f'{stripped_path}: {value}\n')
    print(f"Robustness scores saved to {output_fp}")
    
if __name__ == "__main__":
    main()