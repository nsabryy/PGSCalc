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

def run_analysis(input_fp, index_fp, ks=False, anova=False, levenes=False, bartletts=False, anderson_darling=False):
    output_fp = input_fp.replace('.tsv', '')
    subprocess.run([
        'python3', 'analysis.py',
        '--input', input_fp,
        '--index', index_fp,
        '--output', output_fp,
        '--ks', str(ks),
        '--anova', str(anova),
        '--levenes', str(levenes),
        '--bartletts', str(bartletts),
        '--andar', str(anderson_darling)
    ], check=True)
    
def run_ks_plot(input_fp):
    csv_fp = input_fp.replace('.tsv', '_ks.csv')
    png_fp = input_fp.replace('.tsv', '_ks_graph.png')
    subprocess.run([
        'Rscript', 'ks_plot.r',
        csv_fp,   
        png_fp 
    ], check=True)

def main():
    # USER INPUTS
    base_dir = '/N/project/compgen/PGSCalc/scoring_results/9e09114e_0_0_union/EFO_0004736'
    input_files = find_tsv_files(base_dir)
    index_fp = '/N/project/compgen/shared/resources/1000_genomes/20220422_3202_phased_SNV_INDEL_SV/1000genomes.sequence.index'
    total_files = len(input_files)

    ks = False
    anova = True
    levenes = True
    bartletts = True
    anderson_darling = True

    for i, input_fp in enumerate(input_files, start=1):
        print(f"[{i}/{total_files}] Processing {input_fp}...", end='', flush=True)
        run_analysis(input_fp, index_fp, ks, anova, levenes, bartletts, anderson_darling)
        if ks: run_ks_plot(input_fp)
        print(f"\r[{i}/{total_files}] Processing {input_fp}... Done.")

if __name__ == "__main__":
    main()
