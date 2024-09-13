import argparse
import pandas as pd
from scipy import stats

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
    pop_distributions = {'ALL': []}
    all_scores = []
    with open(scoring_results_fp, 'r') as f:
        # Skip the header line in the scoring file
        header = f.readline()
        for line in f:
            fields = line.strip().split('\t')
            sample_id = fields[1]
            score = float(fields[3])
            
            population = sample_to_population.get(sample_id, "UNKNOWN")
            # if population == "UNKNOWN": 
            #     print(sample_id)
            #     input()
            if population not in pop_distributions:
                pop_distributions[population] = []
            pop_distributions[population].append(score)
            pop_distributions['ALL'].append(score)

    return pop_distributions

def apply_ks_to_population_distributions(population_dict, output_filepath):
    # Prepare a list to store results as dictionaries
    output_filepath = output_filepath + '_ks.csv'

    ks_results = []

    populations = list(population_dict.keys())
    num_populations = len(populations)
    
    for i in range(num_populations):
        for j in range(i + 1, num_populations):
            pop1 = populations[i]
            pop2 = populations[j]
            scores1 = population_dict[pop1]
            scores2 = population_dict[pop2]
            
            # Apply KS test
            ks_statistic, p_value = stats.ks_2samp(scores1, scores2)
            
            # Append results as a dictionary
            ks_results.append({
                "Population Pair": f"{pop1}/{pop2}",
                "KS Statistic": ks_statistic,
                "p-value": p_value
            })
    
    # Convert results to a DataFrame
    df = pd.DataFrame(ks_results)
    
    # Save the DataFrame as a tab-delimited text file
    df.to_csv(output_filepath, sep='\t', index=False)
    
    # print(f"KS test results saved to {output_filepath}")
    return df

def ks_average_vs_population(population_dict, output_filepath):
    output_filepath = output_filepath + '_ks.csv'

    # Prepare a list to store results as dictionaries
    
    ks_results = []

    populations = list(population_dict.keys())
    num_populations = len(populations)
    
    all_scores = population_dict['ALL']
    for i in range(num_populations):
        pop1 = populations[i]
        if pop1 == 'ALL':
            continue

        scores1 = population_dict[pop1]
        # Apply KS test
        ks_statistic, p_value = stats.ks_2samp(scores1, all_scores)
        
        # Append results as a dictionary
        ks_results.append({
            "Population Pair": f"{pop1}",
            "KS Statistic": ks_statistic,
            "p-value": p_value
        })
    
    # Convert results to a DataFrame
    df = pd.DataFrame(ks_results)
    
    # Save the DataFrame as a tab-delimited text file
    df.to_csv(output_filepath, sep='\t', index=False)
    
    # print(f"KS test results saved to {output_filepath}")
    return df

def apply_anova_to_population_distributions(population_dict, output_filepath):
    output_filepath = output_filepath + '_anova.csv'

    # Prepare a list to store the ANOVA results
    anova_results = []

    populations = list(population_dict.keys())
    if 'ALL' in populations:
        populations.remove('ALL')

    # Collect scores for each population
    scores = [population_dict[pop] for pop in populations]
    
    # Apply ANOVA test
    f_statistic, p_value = stats.f_oneway(*scores)
    
    # Append results to a dictionary
    anova_results.append({
        "Test": "ANOVA",
        "F-Statistic": f_statistic,
        "p-value": p_value
    })

    # Convert results to a DataFrame
    df = pd.DataFrame(anova_results)
    
    # Save the DataFrame as a tab-delimited text file
    df.to_csv(output_filepath, sep='\t', index=False)
    
    return df

def apply_levenes_to_population_distributions(population_dict, output_filepath):
    output_filepath = output_filepath + '_levenes.csv'
    
    # Prepare a list to store the Levene's test results
    levenes_results = []
    
    populations = list(population_dict.keys())
    if 'ALL' in populations:
        populations.remove('ALL')

    # Collect scores for each population
    scores = [population_dict[pop] for pop in populations]
    
    # Apply Levene's test
    statistic, p_value = stats.levene(*scores)
    
    # Append results to a dictionary
    levenes_results.append({
        "Test": "Levene's Test",
        "Statistic": statistic,
        "p-value": p_value
    })

    # Convert results to a DataFrame
    df = pd.DataFrame(levenes_results)
    
    # Save the DataFrame as a tab-delimited text file
    df.to_csv(output_filepath, sep='\t', index=False)
    
    return df

def apply_bartletts_to_population_distributions(population_dict, output_filepath):
    output_filepath = output_filepath + '_bartletts.csv'
    
    # Prepare a list to store the Bartlett's test results
    bartletts_results = []
    
    populations = list(population_dict.keys())
    if 'ALL' in populations:
        populations.remove('ALL')

    # Collect scores for each population
    scores = [population_dict[pop] for pop in populations]
    
    # Apply Bartlett's test
    statistic, p_value = stats.bartlett(*scores)
    
    # Append results to a dictionary
    bartletts_results.append({
        "Test": "Bartlett's Test",
        "Statistic": statistic,
        "p-value": p_value
    })

    # Convert results to a DataFrame
    df = pd.DataFrame(bartletts_results)
    
    # Save the DataFrame as a tab-delimited text file
    df.to_csv(output_filepath, sep='\t', index=False)
    
    return df

def apply_anderson_darling_to_population_distributions(population_dict, output_filepath):
    output_filepath = output_filepath + '_anderson_darling.csv'
    
    # Prepare a list to store the Anderson-Darling test results
    ad_results = []
    
    populations = list(population_dict.keys())
    if 'ALL' in populations:
        populations.remove('ALL')

    for pop in populations:
        # Apply Anderson-Darling test
        statistic, critical_values, significance_level = stats.anderson(population_dict[pop])
        
        # Store results
        ad_results.append({
            "Population": pop,
            "Statistic": statistic,
            "Critical Values": critical_values,
            "Significance Level": significance_level
        })

    # Convert results to a DataFrame
    df = pd.DataFrame(ad_results)
    
    # Save the DataFrame as a tab-delimited text file
    df.to_csv(output_filepath, sep='\t', index=False)
    
    return df

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Run KS analysis on population distributions.')
    parser.add_argument('--input', required=True, help='Path to the input TSV file.')
    parser.add_argument('--index', required=True, help='Path to the index file.')
    parser.add_argument('--output', required=True, help='Path to save the results file.')
    parser.add_argument('--ks', required=False, help='Calculate KS scores.')
    parser.add_argument('--anova', required=False, help='Caculate ANOVA scores.')
    parser.add_argument('--levenes', required=False, help='Run Levenes test.')
    parser.add_argument('--bartletts', required=False, help='Run Bartletts test.')
    parser.add_argument('--andar', required=False, help='Run Anderson-Darling test.')

    # Parse arguments
    args = parser.parse_args()
    
    # Read arguments
    scoring_results_fp = args.input
    index_fp = args.index
    output_filepath_data = args.output

    population_dict = get_population_distributions(scoring_results_fp, index_fp)

    if bool(args.ks):
        ks_average_vs_population(population_dict, output_filepath_data)
        # apply_ks_to_population_distributions(population_dict, output_filepath_data)

    if bool(args.anova):
        apply_anova_to_population_distributions(population_dict, output_filepath_data)

    if bool(args.levenes):
        apply_levenes_to_population_distributions(population_dict, output_filepath_data)
        
    if bool(args.bartletts):
        apply_bartletts_to_population_distributions(population_dict, output_filepath_data)
        
    if bool(args.andar):
        apply_anderson_darling_to_population_distributions(population_dict, output_filepath_data)

if __name__ == "__main__":
    main()