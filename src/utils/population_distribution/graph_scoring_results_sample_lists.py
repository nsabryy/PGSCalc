import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import fnmatch
def generate_graph_title(file_path):
    # Extract file name from the file path
    file_name = os.path.basename(file_path)
    
    # Check if "union" is in the file name before the extension
    if 'union' in file_name:
        snp_type = "common snps"
    else:
        snp_type = "all snps"
    
    # Extract the PGS ID (first part of the file name)
    pgs_id = file_name.split('_')[0]
    
    # Construct the graph title
    title = f"{pgs_id} Distribution by Cohort, {snp_type}"
    
    return title

def graph_by_cohort(samplelist1, samplelist2, scoring_results_fp):
    # Read cohort sample lists
    with open(samplelist1[1], 'r') as file:
        samples1 = file.read().splitlines()
    with open(samplelist2[1], 'r') as file:
        samples2 = file.read().splitlines()

    # Format the title
    title = generate_graph_title(scoring_results_fp)
        
    # Map samples to cohort
    s1_scores = []
    s2_scores = []

    # Read the scoring results file line by line
    fp = scoring_results_fp
    with open(fp, 'r') as f:
        # Skip the header
        header = f.readline()
        for line in f:
            fields = line.strip().split('\t')
            sample = fields[1]
            score = float(fields[3])
            if sample in samples1:
                s1_scores.append(score)
            elif sample in samples2:
                s2_scores.append(score)

    # Convert lists to numpy arrays for consistency
    s1_scores = np.array(s1_scores)
    s2_scores = np.array(s2_scores)

    try:
        bins = np.linspace(min(min(s1_scores), min(s2_scores)), max(max(s1_scores), max(s2_scores)), 30)
    except Exception as e: 
        print(e)
        return

    # Create histograms
    hist_s1, bins_s1 = np.histogram(s1_scores, bins=bins)
    hist_s2, bins_s2 = np.histogram(s2_scores, bins=bins)

    # Plot histograms
    plt.hist(bins_s1[:-1], bins_s1, weights=hist_s1, alpha=0.5, label=samplelist1[0], color='blue')
    plt.hist(bins_s2[:-1], bins_s2, weights=hist_s2, alpha=0.5, label=samplelist2[0], color='red')

    # Add labels and legend
    plt.xlabel('Score')
    plt.ylabel('Density')
    plt.title(title)    
    plt.legend(loc='upper right')

    # Show the plot
    # plt.show()
    output_fp = fp.replace(".tsv", "_cohort_distribution.png")
    plt.savefig(output_fp)

def main():
    samplelist1 = ['EUR', '/N/project/compgen/PGSCalc/input_data/cohort_sample_lists/eur_samples.txt']
    samplelist2 = ['AFR', '/N/project/compgen/PGSCalc/input_data/cohort_sample_lists/afr_samples.txt']
    # for root, dirs, files in os.walk(root_dir):
    #     dirs[:] = [d for d in dirs if d not in ignore_dirs]
    #     for file in fnmatch.filter(files, '*.tsv'):
    #         file_path = os.path.join(root, file)
    #         print(f"generating cohort specific graph for {file_path}")
    #         graph_by_cohort(samplelist1, samplelist2, file_path)
    scoring_results_fp = '/N/project/compgen/PGSCalc/scoring_results/f990859b_0_0_union/PGS000001/PGS000001_f990859b_0_0_union.tsv'
    graph_by_cohort(samplelist1, samplelist2, scoring_results_fp)

if __name__ == "__main__":
    main()
