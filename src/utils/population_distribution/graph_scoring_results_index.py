import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
    pop_distributions = {}
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

    return pop_distributions

def graph_population_histograms(scoring_results_fp, index_fp, line_density =True):
    individual_histograms = get_population_distributions(scoring_results_fp, index_fp)
    # Choose whether to plot histograms or line density plots
    if not line_density:
        for pop, scores in individual_histograms.items():
            plt.hist(scores, bins=30, alpha=0.5, label=f'{pop}', density=False)
        ylabel = 'Frequency'
        output_fp = scoring_results_fp.replace(".tsv", "_population_histogram.png")
    else:
        for pop, scores in individual_histograms.items():
            sns.kdeplot(scores, label=f'{pop}', fill=False)  # Line density plot
        ylabel = 'Density'
        output_fp = scoring_results_fp.replace(".tsv", "_population_density_plot.png")
    
    plt.xlabel('Score')
    plt.ylabel(ylabel)
    plt.title('Population-specific Score Distribution')
    plt.legend()
    
    # Save the plot
    print(output_fp)
    plt.savefig(output_fp)
    plt.close()

def main():
    scoring_results_fp = '/N/project/compgen/PGSCalc/scoring_results/9e09114e_0_0_union/EFO_0004736/PGS004722/PGS004722_9e09114e_0_0_union.tsv'
    index_fp = '/N/project/compgen/shared/resources/1000_genomes/20220422_3202_phased_SNV_INDEL_SV/1000genomes.sequence.index'
    graph_population_histograms(scoring_results_fp, index_fp)

if __name__ == "__main__":
    main()