import re
import hashlib
import os
import logging
import pandas as pd
import pickle
import json
import argparse
from mpi4py import MPI
from VCFQuery import VCFQuery
from PGScore import PGScore 
import traceback

def parse_config():
    """Parse the configuration file argument."""
    parser = argparse.ArgumentParser(description="Download PGS scoring files and query sample data.")
    parser.add_argument('-c', '--config', required=True, help="Path to the configuration file.")
    parser.add_argument('--job_id', required=True, help="Job ID for this run.")
    return parser.parse_args()

def load_config(config_path):
    """Load the configuration file."""
    with open(config_path, 'r') as f:
        config = json.load(f)
    if type(config['sample_ids'])==str and os.path.isfile(config['sample_ids']):
        with open(config['sample_ids'], 'r') as sample_file:
            config['sample_ids']=sample_file.read().splitlines()
            print("Parsing file provided for sampleIDs")
    return config

def setup_logging(job_id):
    """Set up logging with the given job ID."""
    log_dir = f'logs/{job_id}'
    os.makedirs(log_dir, exist_ok=True)
    log_filename = f'{log_dir}/pgscalc_{job_id}.log'
    logging.basicConfig(filename=log_filename, level=logging.INFO, 
                        format='%(asctime)s %(levelname)s: %(message)s')
    # logging.info(f"Logging setup complete. Log file: {log_filename}")
    return log_filename

def generate_output_filename(input_dir, sample_ids):
    """Generate a unique file name based on input directory and optional input parameters: samples, imputed, filtered"""
    try:
        input_dir_normalized = re.sub(r'[\\/:*?"<>|]', '_', input_dir.strip('/').replace(' ', '_'))
        
        if sample_ids:
            # Generate a hash from the sample IDs to ensure uniqueness
            sample_ids = ','.join(sorted(sample_ids))
            sample_ids_hash = hashlib.md5(sample_ids.encode()).hexdigest()[:8]
            logging.info(f"Unique subset ID for samples provided is: {sample_ids_hash}")
            # Append the hash to the normalized input directory name
            input_dir_normalized = f"{input_dir_normalized}_subset_{sample_ids_hash}"

        return input_dir_normalized
    except Exception as e:
        logging.error(f"Error generating output filename: {e}")
        raise

def format_call(call, ref, alts):
    """Change numeric GT to named GT."""
    try:
        if call is None:
            logging.warning("Genotype call is missing.")
            return None

        if ref is None:
            logging.warning("Reference allele is missing.")
            return None

        alts = alts if alts is not None else []
        allele_list = [ref] + list(alts)

        formatted_call = tuple(
            allele_list[c] if c is not None and c != '.' else c for c in call
        )
        return formatted_call
    except Exception as e:
        logging.error(f"Error in format_call: {e}")
        raise

def get_combined_scores(config):
    """Pull combined scores dictionary from PGScore."""
    try:
        if 'scoring_file_directory' in config['PGScore'] and config['PGScore']['scoring_file_directory']:
            logging.info(f"Loading local PGS files: {config['PGScore']['scoring_file_directory']}")
            score_package = PGScore.ScorePackage(
                scoring_file_directory=config['PGScore']['scoring_file_directory']
            )
        else:
            logging.info(f"Starting download for PGS IDs: {config['PGScore']['pgs_ids']}")
            score_package = PGScore.ScorePackage(
                pgs_ids=config['PGScore']['pgs_ids'],  
                genome_build=config['PGScore']['genome_build'],
                override=config['PGScore']['override']
            )
        return score_package
    except Exception as e:
        logging.error(f"Error in configuration or downloading scores: {e}")
        return None

def process_records(ref, chrom=[], samples=None):
    """Process SNPs and dynamically build the results dictionary."""
    x = 1
    logging.info("Starting processing of records")
    results_dict = {}
    good_alleles = {"A", "C", "T", "G"}
    # Query the SNPs
    for record in ref.query_imputed_snps(
        chrom=chrom, 
        samples=samples
    ):
        chrom = record.chrom
        pos = record.pos

        for sample in record.samples:
            if samples and sample not in samples:
                continue

            key = (chrom, pos)

            if key not in results_dict:
                results_dict[key] = {
                    'A': 0,
                    'C': 0,
                    'T': 0,
                    'G': 0,
                    'unknown': 0
                }
            try:
                call = format_call(record.samples[sample]['GT'], record.ref, record.alts)
            except KeyError as e:
                logging.error(f"Call missing for sample {sample} at {chrom}:{pos}")
                continue

            if call is None:
                continue
                
            for allele in call:
                if allele in good_alleles:
                    results_dict[key][allele] += 1
                else:
                    results_dict[key]['unknown'] += 1
            
            logging.debug(f"Updated results for {key}: {results_dict[key]}")

            if x == 1:
                print(results_dict[key])
                x += 1
    return results_dict  

def sequential_main(config,args):
    setup_logging(args.job_id)
    logging.info("Job started")
    input_dir = config['input_dir']
    ref = VCFQuery.VCFQuery(vcf_location=input_dir)
    logging.info("VCFQuery object made")
    sample_ids = config['sample_ids']
    if not sample_ids:
        sample_ids = ref._samples

    results_dict = process_records(ref, config['chrom'], samples=sample_ids)
    
    results_df = pd.DataFrame.from_dict(results_dict, orient='index')
    results_df.to_csv("/N/project/compgen/PGSCalcimputation_analysis_results/testing.csv", sep='\t', index=False)

def main():
    args = parse_config()
    config = load_config(args.config)
    if config["parallel"]:
        try:
            multiprocessing_main(config,args)
        except ValueError as e:
            logging.error(e)
            sequential_main(config,args)
    else: 
        sequential_main(config,args)

if __name__ == "__main__":
    print("Starting job!")
    main()


# def multiprocessing_main(config,args):
#     try:
#         comm = MPI.COMM_WORLD
#         rank = comm.Get_rank()
#         size = comm.Get_size()

#         if rank == 0:
#             input_dir = config['input_dir']
#             ref = VCFQuerya(vcf_location=input_dir)
#             ref_pkl = pickle.dumps(ref)
            
#             sample_ids_fp = config['sample_ids']
#             with open(sample_ids_fp, 'r') as fp:
#                 sample_ids = fp.read().splitlines()
#             if not sample_ids:
#                 sample_ids = ref._samples

#             # # Create PGScore object to fetch specific SNPs
#             # score_package = get_combined_scores(config)   # Get SNPs related to this PGS ID

#             vcf_sources = ref.query_sources(chrom=config['chrom'], samples=sample_ids)
#             logging.info(f"VCF sources: {vcf_sources}")
#             # Distribute the VCF sources among the processes
#             vcf_sources_split = [vcf_sources[i::size] for i in range(size)]
#         else:
#             ref_pkl = None
#             sample_ids = None
#             vcf_sources_split = None

#         # Broadcast the reference object and sample IDs to all processes
#         ref_pkl = comm.bcast(ref_pkl, root=0)
#         sample_ids = comm.bcast(sample_ids, root=0)
#         vcf_sources_split = comm.scatter(vcf_sources_split, root=0)

#         # Deserialize the reference object
#         ref = pickle.loads(ref_pkl)

#         ##### CALC #####
#         results_df = evaluate_imputed(ref, vcf_sources_split, sample_ids)

#         # Gather the results from all processes
#         all_results = comm.gather(results_df, root=0)

#         if rank == 0:
#             # Combine the results into a single DataFrame
#             combined_results_df = pd.concat(all_results)

#             # Save the results to a file
#             output_dir = config['output_dir']
#             os.makedirs(output_dir, exist_ok=True)
#             output_filename = generate_output_filename(input_dir, sample_ids)
#             output_file = os.path.join(output_dir, f"{output_filename}.tsv")
#             combined_results_df.to_csv(output_file, sep='\t', index=False)
#             logging.info(f"Results saved to {output_file}")

#     except Exception as e:
#         logging.error(f"Error in main: {e}")
#         logging.error(traceback.format_exc())
