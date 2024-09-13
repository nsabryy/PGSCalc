import pandas as pd
import ast

def calculate_entropy(snp_calls):
    """
    Returns a dictionary with chrom:pos as keys and call distributions as values.
    """
    entropy_dict = {}

    for chrom, pos, ref, alt, call in snp_calls:
        key = f"{chrom}:{pos}"
        if key not in entropy_dict:
            entropy_dict[key] = {
                'ref': ref,
                'alt': alt,
                'homozygous_ref': 0,
                'homozygous_alt': 0,
                'heterozygous': 0
            }

        if call == ref + ref:
            entropy_dict[key]['homozygous_ref'] += 1
        elif call == alt + alt:
            entropy_dict[key]['homozygous_alt'] += 1
        else:
            entropy_dict[key]['heterozygous'] += 1

    return entropy_dict


def update_sample_data(df):
    """
    Updates the DataFrame to reflect the SNP calls per sample and calculates the entropy distributions.
    """
    overall_entropy = {}

    for index, row in df.iterrows():
        snp_calls = parse_snp_column(row['SNP'])
        sample_entropy = calculate_entropy(snp_calls)
        
        # Update the overall entropy distribution
        for key, value in sample_entropy.items():
            if key not in overall_entropy:
                overall_entropy[key] = {
                    'ref': value['ref'],
                    'alt': value['alt'],
                    'homozygous_ref': 0,
                    'homozygous_alt': 0,
                    'heterozygous': 0
                }
            else:
                # Ensure the ref/alt alleles match before combining counts
                if overall_entropy[key]['ref'] != value['ref'] or overall_entropy[key]['alt'] != value['alt']:
                    # If ref/alt are swapped, switch counts
                    overall_entropy[key]['homozygous_ref'] += value['homozygous_alt']
                    overall_entropy[key]['homozygous_alt'] += value['homozygous_ref']
                else:
                    overall_entropy[key]['homozygous_ref'] += value['homozygous_ref']
                    overall_entropy[key]['homozygous_alt'] += value['homozygous_alt']

            overall_entropy[key]['heterozygous'] += value['heterozygous']

    return overall_entropy


def parse_snp_column(snp_string):
    """
    Parses the complex SNP column string to extract each chrom:pos and its alleles.
    Returns a list of tuples in the format: [(chrom, pos, ref, alt, call), ...]
    """
    snp_data = ast.literal_eval(snp_string)  # Safely evaluates the string as a Python object (list of lists/tuples)
    parsed_data = []

    for snp_set in snp_data:
        if isinstance(snp_set, list):
            for chrom, pos, ref, alt_list, call in snp_set:
                # Handle cases where there are multiple alternative alleles
                if len(alt_list) > 1:
                    print(f"Multiple alts encountered at {chrom}:{pos} - skipping this SNP.")
                    continue
                alt = alt_list[0]
                parsed_data.append((chrom, pos, ref, alt, ''.join(call)))
        else:
            chrom, pos, ref, alt_list, call = snp_set
            if len(alt_list) > 1:
                print(f"Multiple alts encountered at {chrom}:{pos} - skipping this SNP.")
                continue
            alt = alt_list[0]
            parsed_data.append((chrom, pos, ref, alt, ''.join(call)))

    return parsed_data


def main():
    # Load the TSV file into a DataFrame
    df = pd.read_csv('/N/project/compgen/PGSCalc/scoring_results/testing/N_project_biobank_imputed-from-dongbing_reheader_subset_f3c49daf_0_0_union/PGS000001/PGS000001_N_project_biobank_imputed-from-dongbing_reheader_subset_f3c49daf_0_0_union.tsv', sep='\t')
    
    # Calculate and update the entropy distribution
    entropy_results = update_sample_data(df)
    entropy_df = pd.DataFrame(entropy_results).T 
    # Save results to a new TSV file
    entropy_df.to_csv('/N/project/compgen/PGSCalc/scoring_results/testing/N_project_biobank_imputed-from-dongbing_reheader_subset_f3c49daf_0_0_union/PGS000001/PGS000001_N_project_biobank_imputed-from-dongbing_reheader_subset_f3c49daf_0_0_union_imputed_analysis.tsv', sep='\t')

if __name__ == "__main__":
    main()
