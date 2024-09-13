from mpi4py import MPI
import argparse
import logging
import hashlib
import os
import re
import json
import pandas as pd
import pickle
import numpy as np
import matplotlib.pyplot as plt
import heapq
from scipy.stats import zscore
from sklearn.mixture import GaussianMixture
import sqlite3
import sys
from VCFQuery import VCFQuery
from PGScore import PGScore

# @contextmanager
def suppress_logging(level=logging.CRITICAL):
    logger = logging.getLogger()
    previous_level = logger.level
    logger.setLevel(level)
    try:
        yield
    finally:
        logger.setLevel(previous_level)

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
    if type(config['VCFQuery']['sample_ids'])==str and os.path.isfile(config['VCFQuery']['sample_ids']):
        with open(config['VCFQuery']['sample_ids'], 'r') as sample_file:
            config['VCFQuery']['sample_ids']=sample_file.read().splitlines()
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

def generate_output_filename(config):
    """Generate a unique file name based on input directory and optional input parameters: samples, imputed, filtered"""
    # Combine all input directories into a single string and generate a hash
    input_dirs = ','.join(sorted(config['VCFQuery']['input_dir']))
    input_dir_hash = hashlib.md5(input_dirs.encode()).hexdigest()[:8]
    input_dir_normalized = f"{input_dir_hash}"

    if config['VCFQuery']['sample_ids']:
        # Generate a hash from the sample IDs to ensure uniqueness
        sample_ids = ','.join(sorted(config['VCFQuery']['sample_ids']))
        sample_ids_hash = hashlib.md5(sample_ids.encode()).hexdigest()[:8]
        logging.info(f"Unique subset ID for samples provided is: {sample_ids_hash}")
        # Append the hash to the normalized input directory name
        input_dir_normalized = f"{input_dir_normalized}_subset_{sample_ids_hash}"

    # Get the imputed and filtered flags from the config
    imputed_flag = 1 if config['VCFQuery'].get('remove_imputed', False) else 0
    filtered_flag = 1 if config['VCFQuery'].get('remove_filtered', False) else 0

    # Append the flags to the normalized input directory name
    output_filename = f"{input_dir_normalized}_{imputed_flag}_{filtered_flag}"
    if config['PGSCalc'].get('calculate_overlap_only', False): 
        output_filename = output_filename + "_union"
    return output_filename

def calculate_snp_score(call, effect_allele, effect_weight):
    """Calculate the SNP score based on the call and effect allele."""
    if call[0] not in {'A', 'C', 'T', 'G'} and call[1] not in {'A', 'C', 'T', 'G'}:
        # logging.error(f"Call not processed; call: {call}")
        return 0
    return sum(effect_weight for allele in call if allele == effect_allele)

def format_call(call, ref, alts):
    """Change numeric GT to named GT."""
    # Check if genotype call is missing
    if call is None:
        logging.error("Genotype call is missing.")
        return None

    # Check if reference allele is missing
    if ref is None:
        logging.error("Reference allele is missing.")
        return None

    # Handle alternate alleles if they are missing
    alts = alts if alts is not None else []

    # Create a list of alleles
    allele_list = [ref] + list(alts)

    # Format the call
    formatted_call = tuple(
        allele_list[c] if c is not None and c != '.' else c for c in call
    )
    return formatted_call

def process_snps(combined_scores, ref, score_metadata, samples=None, remove_imputed=False, remove_filtered=False, include_snps=False, only_imputed=False, remove_at_gc=False):
    """Process SNPs and dynamically build the results dictionary."""
    results_dict = {}

    # Extract the SNPs to query in the format {chrom: [pos]}
    snps_to_query = {chrom: [int(key) for key in pos_dict.keys()] for chrom, pos_dict in combined_scores.items()}

    all_samples = set(ref._samples.keys())
    last_snp = None
    # Query the SNPs
    for query_result in ref.query_snps(
        snps=snps_to_query, 
        samples=samples, 
        remove_imputed=remove_imputed, 
        remove_filtered=remove_filtered
    ):
        filename, record = query_result
        chrom = record.chrom
        if 'chr' in chrom:
            chrom = chrom[3:]
        pos = record.pos
        ref_set = set(record.ref)
        alt_set = set(record.alts) if record.alts != "<NON_REF>" else set('.')

        if not last_snp or last_snp!=(chrom,pos):   
            last_snp = (chrom, pos)
            sample_check_per_score = {score['score']: all_samples.copy() for score in combined_scores[chrom][pos]}

        # logging.info(f"Received record for {chrom}:{pos}, we have {len(record.samples)} samples")
        if remove_at_gc:
            if ('A' in ref_set and 'T' in alt_set) or ('T' in ref_set and 'A' in alt_set):
                logging.info(f"A-T SNP found for snp {chrom}:{pos}. Ignoring for scoring.")
                continue
            elif ('G' in ref_set and 'C' in alt_set) or ('C' in ref_set and 'G' in alt_set):
                logging.info(f"G-C SNP found for snp {chrom}:{pos}. Ignoring for scoring.")
                continue

        if chrom not in combined_scores or pos not in combined_scores[chrom]:
            continue

        for scoring_info in combined_scores[chrom][pos]:
            pgs_id = scoring_info['score']

            # Allele verification
            eff_set = set(scoring_info['effect_allele'])
            oth_set = set(scoring_info['other_allele']) if scoring_info['other_allele'] else set('.')

            if not allele_verification(ref_set, alt_set, eff_set, oth_set):
                logging.warning(f"In file {filename} for score {pgs_id} Input reference or ALT allele mismatch with VCF at {chrom}:{pos}. SCOREFILE: Effect allele: {eff_set}, Other Allele: {oth_set}; VCF FILE: Ref: {ref_set}, Alt: {alt_set}")
                continue

            for sample in record.samples:
                if samples and sample not in samples:
                    continue
                
                # Verify we have only seen this sample 1x for each score
                if sample not in sample_check_per_score[pgs_id]:
                    # logging.warning(f"Received repeat record info for sample {sample} for {chrom}:{pos}.")
                    continue

                sample_check_per_score[pgs_id].remove(sample)

                key = (pgs_id, sample)

                if key not in results_dict:
                    results_dict[key] = {
                        'PGS_ID': pgs_id,
                        'SAMPLE': sample,
                        'TOTAL_ALLELES_QUERIED': score_metadata.get(pgs_id, 0),  
                        'SCORE': 0,
                        'ALLELES_FOUND': 0,
                        'NUMBER_IMPUTED': 0
                    }
                    if include_snps:
                        results_dict[key]['SNP'] = []
                try:
                    call = format_call(record.samples[sample]['GT'], record.ref, record.alts)
                except KeyError as e:
                    logging.error(f"Call missing for sample {sample} at {chrom}:{pos}")
                    continue

                if call is None:
                    continue

                sample_row = results_dict[key]

                # Update alleles found for sample
                sample_row['ALLELES_FOUND'] += 1
                # Add the snp id if requested
                if include_snps:
                    sample_row['SNP'].append((record.chrom, record.pos, record.ref, record.alts, call))
                # Update number of imputed
                if 'IMPUTED' in record.info and record.info['IMPUTED']:
                    sample_row['NUMBER_IMPUTED'] += 1
                # Update score
                sample_row['SCORE'] += calculate_snp_score(
                    call=call, 
                    effect_allele=scoring_info['effect_allele'], 
                    effect_weight=scoring_info['effect_weight']
                    )

    return results_dict   

def allele_verification(ref, alt, eff, oth):
    """Verify alleles between the VCF and the scoring file."""
    if "<NON_REF>" in alt:
        if oth == set(['.']):
            if not ref == eff:
                print("Not enough information ... both other_allele and alternative allele information missing and effect doesn't match other.")
                return False
        else:
            if not (eff == ref or oth == ref):
                return False
    else:
        if oth == set(['.']):
            if not ((eff == ref) or (eff == alt)):
                return False
        else:
            if not ((ref == eff and alt == oth) or (ref == oth and alt == eff)):
                return False
    return True

def calculate_and_save_summary_statistic(results_df, config, output_filename, pgs_id):
    """Calculate and save summary statistics for a single PGS ID."""
    # Calculate % match
    total_alleles_queried = results_df['TOTAL_ALLELES_QUERIED'].sum()
    alleles_found = results_df['ALLELES_FOUND'].sum()
    match_percentage = (alleles_found / total_alleles_queried) * 100

    # Get top 10% scores
    top_10_percent_threshold = results_df['SCORE'].quantile(0.9)
    top_10_percent_samples = results_df[results_df['SCORE'] >= top_10_percent_threshold]

    # Create directories if they do not exist
    summary_dir = os.path.join(config['PGSCalc']['output_dir'], output_filename, pgs_id)
    os.makedirs(summary_dir, exist_ok=True)

    # Save summary to a text file
    summary_file = os.path.join(summary_dir, f"{pgs_id}_{output_filename}_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Configuration Details:\n")
        for section, section_config in config.items():
            f.write(f"[{section}]\n")
            for key, value in section_config.items():
                f.write(f"{key}: {value}\n")
            f.write("\n")

        f.write(f"Summary Statistics for PGS ID: {pgs_id}\n")
        f.write(f"Total Alleles Queried: {total_alleles_queried}\n")
        f.write(f"Alleles Found: {alleles_found}\n")
        f.write(f"Match Percentage: {match_percentage:.2f}%\n")
        f.write(f"Top 10% Score Threshold: {top_10_percent_threshold}\n")
        f.write(f"Top 10% Samples:\n")
        f.write(top_10_percent_samples[['SAMPLE', 'SCORE']].to_string(index=False))

    # Calculate z-scores and add to DataFrame
    results_df.loc[:, 'RISK_SCORE'] = zscore(results_df['SCORE'])

    # Save the distribution type to the summary file
    with open(summary_file, 'a') as f:
        f.write(f"\nDistribution Type: unimodal\n")

    # Generate and save a distribution plot of scores
    plot_dir = os.path.join(config['PGSCalc']['output_dir'], output_filename, pgs_id)
    os.makedirs(plot_dir, exist_ok=True)
    
    plt.figure(figsize=(10, 6))
    plt.hist(results_df['SCORE'], bins=30, density=True, alpha=0.75)
    plt.title(f"Distribution of Scores for PGS ID: {pgs_id}")
    plt.xlabel("Score")
    plt.ylabel("Density")
    plt.grid(True)
    plt.savefig(os.path.join(plot_dir, f"{pgs_id}_{output_filename}_score_distribution.png"))
    plt.close()
    
    logging.info(f"Summary statistics saved to: {summary_file}")

def save_results(results_df, config, output_filename, pgs_id):
    """Convert results dictionary to DataFrame and save to a file."""    
    output_dir = os.path.join(config['PGSCalc']['output_dir'], output_filename, pgs_id)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{pgs_id}_{output_filename}.tsv")
    results_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Saved scoring output to: '{output_file}'")

def get_combined_scores(config):
    """Pull combined scores dictionary from PGScore."""
    try:
        if 'scoring_file_directory' in config['PGScore'] and config['PGScore']['scoring_file_directory']:
            logging.info(f"Loading local PGS files: {config['PGScore']['scoring_file_directory']}")
            score_package = PGScore.ScorePackage(
                scoring_file_directory=config['PGScore']['scoring_file_directory']
            )
        else:
            logging.info(f"Starting download for PGS.")
            score_package = PGScore.ScorePackage(
                pgs_ids=config['PGScore'].get('pgs_ids', None), 
                trait_id=config['PGScore'].get('trait_id', None),
                trait_include_children=config['PGScore'].get('trait_include_children', 1), 
                genome_build=config['PGScore']['genome_build'],
                override=config['PGScore']['override']
            )
        return score_package

    except Exception as e:
        logging.error(f"Error in configuration or downloading scores: {e}")
        return None

# TODO we don't need this Verify data integrity
def verify_data_integrity(input_dict, output_dicts):
    """
    Verifies that all data in the input_dict is still present in the output_dicts.

    Args:
    input_dict (dict): The original input dictionary structured as {chrom: {pos: [info]}}.
    output_dicts (list): A list of output dictionaries, each structured similarly to the input_dict.

    Returns:
    bool: True if all data is intact, False otherwise.
    """
    for chrom, positions in input_dict.items():
        for pos, info in positions.items():
            # Check if this chrom and pos exists in all output dicts
            found = False
            for output_dict in output_dicts:
                if chrom in output_dict and pos in output_dict[chrom]:
                    # Verify the data is identical
                    if output_dict[chrom][pos] == info:
                        found = True
                        break
            
            if not found:
                logging.warning(f"Missing or altered data at chrom {chrom}, pos {pos}")
                return False

    return True

def balance_distribution(number_of_processes, entire_dict, buffer=100):
    """Balances tasks across available processes."""
    # Step 1: Calculate the total size of all key-value pairs
    total_size = sum(len(list(sub_dict.values())[0]) for sub_dict in entire_dict.values() if sub_dict.values())
    if number_of_processes == 0:
        number_of_processes = 1
    avg_load = total_size // number_of_processes

    # Step 2: Sort the keys by the size of their sub-dictionaries in descending order
    sorted_items = sorted(entire_dict.items(), key=lambda item: len(list(item[1].values())[0]), reverse=True)

    # Step 3: Initialize process loads and result dictionaries
    process_loads = [0] * number_of_processes
    result_dicts = [{} for _ in range(number_of_processes)]

    # Use a min-heap to keep track of the process loads
    heap = [(0, i) for i in range(number_of_processes)]
    heapq.heapify(heap)

    # Step 4: Distribute the sub-dictionaries
    for key, sub_dict in sorted_items:
        sub_size = len(list(sub_dict.values())[0])

        # Find the process with the minimum load
        current_load, process_index = heapq.heappop(heap)

        # If adding the entire sub-dictionary exceeds the average load + buffer, split it
        if current_load + sub_size > avg_load + buffer:
            split_size = (current_load + sub_size) - avg_load
            split_sub_dict = {list(sub_dict.keys())[0]: list(sub_dict.values())[0][:split_size]}
            remaining_sub_dict = {list(sub_dict.keys())[0]: list(sub_dict.values())[0][split_size:]}

            # Add the split sub-dictionary to the current process
            result_dicts[process_index][key] = split_sub_dict
            process_loads[process_index] += len(split_sub_dict[list(split_sub_dict.keys())[0]])
            heapq.heappush(heap, (process_loads[process_index], process_index))

            # Find the next process with the minimum load
            current_load, process_index = heapq.heappop(heap)
            result_dicts[process_index][key] = remaining_sub_dict
            process_loads[process_index] += len(remaining_sub_dict[list(remaining_sub_dict.keys())[0]])
        else:
            # Add the entire sub-dictionary to the current process
            result_dicts[process_index][key] = sub_dict
            process_loads[process_index] += sub_size

        # Push the updated load back into the heap
        heapq.heappush(heap, (process_loads[process_index], process_index))

    # Step 5: Resize if necessary
    # Redistribute if the difference in loads exceeds the buffer
    max_load = max(process_loads)
    min_load = min(process_loads)

    if max_load - min_load > buffer:
        # Flatten all sub-dictionaries into a single list of items
        all_items = [(k, v) for subdict in result_dicts for k, v in subdict.items()]
        all_items = sorted(all_items, key=lambda item: len(list(item[1].values())[0]), reverse=True)

        # Reinitialize the process loads and result dictionaries
        process_loads = [0] * number_of_processes
        result_dicts = [{} for _ in range(number_of_processes)]
        heap = [(0, i) for i in range(number_of_processes)]
        heapq.heapify(heap)

        # Redistribute the sub-dictionaries
        for key, sub_dict in all_items:
            sub_size = len(list(sub_dict.values())[0])

            # Find the process with the minimum load
            current_load, process_index = heapq.heappop(heap)

            # Add the entire sub-dictionary to the current process
            result_dicts[process_index][key] = sub_dict
            process_loads[process_index] += sub_size

            # Push the updated load back into the heap
            heapq.heappush(heap, (process_loads[process_index], process_index))
    # Run the integrity check
    # if not verify_data_integrity(entire_dict, result_dicts):
    #     logging.warning("Warning: Data integrity check failed!")
    # else:
    #     logging.info("Data integrity check passed.")

    return result_dicts

def get_pkl_score_distribution(size, combined_scores):
    distributed_scores_pkl = []
    num_partitions = size-1
    size_check = True
    max_allowed_size = 2**31 - 1

    while size_check:
        distributed_scores = balance_distribution(number_of_processes=num_partitions, entire_dict=combined_scores)
        size_check = False 
        distributed_scores_pkl = []

        for partition in distributed_scores:
            part_pkl = pickle.dumps(partition)
            if len(part_pkl) > max_allowed_size:
                logging.warning("Scores to distribute too large. Attempting resize.")
                num_partitions *= 2  
                size_check = True  
                break 
            else:
                distributed_scores_pkl.append(part_pkl)

    for i in range(0, len(distributed_scores_pkl), size):
        distributed_scores_pkl.insert(i, None)

    return distributed_scores_pkl

def process_and_save_score(score, ref, config, output_filename):
    """Process and save the score for a single PGS ID dynamically."""
    results_dict = process_snps(score, ref, config['VCFQuery']['sample_ids'], remove_imputed=config['VCFQuery']['remove_imputed'], remove_filtered=config['VCFQuery']['remove_filtered'])
    results_df = pd.DataFrame.from_dict(results_dict, orient='index')
    calculate_and_save_summary_statistic(results_df, config, output_filename, score.pgs_id)
    save_results(results_df, config, output_filename, score.pgs_id)

# TODO add overlap flags 
def multiprocessing_main(config, args):
    """Multithreaded PGSCalc pipeline."""
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if size == 1:
        raise ValueError("Number of ranks is 1. Defaulting to sequential calculation.")

    setup_logging(args.job_id)
    # logging.info(f"Rank {rank} started")

    # STEP 1: Distribute data
    if rank == 0:
        logging.info(f"Starting computation. Using {size} ranks.")
        logging.info(f"Rank {rank} initializing VCFQuery and PGScore objects.")

        # Pull scoring schema
        score_package = get_combined_scores(config)
        num_of_snps = sum(len(list(sub_dict.values())[0]) for sub_dict in score_package.combined_scores.values() if sub_dict.values())

        if num_of_snps == 0:
            logging.error("No SNPs to query. Check if accessions are valid, or attempt PGScore override if files are corrupted.")
            logging.error("Aborting all ranks.")
            comm.Abort()

        logging.info(f"Beginning query for {sum(len(list(sub_dict.values())) for sub_dict in score_package.combined_scores.values() if sub_dict.values())} number of SNPs")

        # ensure everything is processed & serialized
        ref = VCFQuery.VCFQuery(vcf_locations=config['VCFQuery']['input_dir'], override=config['VCFQuery']['override'])

        # TODO too big when ref big
        # ref_pkl = pickle.dumps(ref)

        # Distribute scores
        logging.info("Balancing scores.")
        distributed_scores_pkl = get_pkl_score_distribution(size=size, combined_scores=score_package.combined_scores)

        # Generate score metadata
        score_metadata = {score.pgs_id: sum(len(positions) for positions in score.snps.values()) for score in score_package.scores}
        score_metadata_pkl = pickle.dumps(score_metadata)

    else:
        distributed_scores_pkl = None
        score_metadata_pkl = None

    # Broadcast and scatter the data
    score_metadata_pkl = comm.bcast(score_metadata_pkl, root=0)
    scattered_scores_pkl = comm.scatter(distributed_scores_pkl, root=0)

    # while True:
    #     try:
    #         while distributed_scores_pkl:
    #             if rank == 0:
    #                 chunk_to_scatter = distributed_scores_pkl[:size]
    #                 distributed_scores_pkl = distributed_scores_pkl[size:]

    #                 if len(chunk_to_scatter) < size:
    #                     chunk_to_scatter += [None] * (size - len(chunk_to_scatter))
    #             else:
    #                 chunk_to_scatter = None

    #             scattered_scores_pkl = comm.scatter(chunk_to_scatter, root=0)

    #         break  

    #     except OverflowError:
    #         if rank == 0:
    #             next_size *= 2
    #             logging.warning(f"Overflow error encountered during scatter. Increasing the number of partitions to {next_size}.")
    #             distributed_scores_pkl = get_pkl_score_distribution(size=next_size, combined_scores=score_package.combined_scores)

    #         distributed_scores_pkl = comm.bcast(distributed_scores_pkl, root=0)


    # STEP 2: Process data
    score_metadata = pickle.loads(score_metadata_pkl)
    chrom_chunk = pickle.loads(scattered_scores_pkl) if scattered_scores_pkl else {}
    with suppress_logging():
        ref = VCFQuery.VCFQuery(vcf_locations=config['VCFQuery']['input_dir'], override=config['VCFQuery']['override'])
    if rank == 0:
        try:
            # Gather results from other ranks
            results = comm.gather(None, root=0)
            logging.info("Gathered all results.")
        except Exception as e:
            logging.error(f"Error during gathering results: {e}")


    else:
        if not chrom_chunk:
            results_dict = {}
        else:
            if config['PGSCalc']['calculate_overlap_only']:
                chrom_chunk = ref.snps_for_cohort(chrom_chunk, samples=config['VCFQuery']['sample_ids'])
            results_dict = process_snps(chrom_chunk, ref, score_metadata,
                                        samples=config['VCFQuery']['sample_ids'],
                                        remove_imputed=config['VCFQuery']['remove_imputed'],
                                        remove_filtered=config['VCFQuery']['remove_filtered'],
                                        include_snps=config['PGSCalc']['incude_snps_per_sample']
                                        )
        results_dict_pkl = pickle.dumps(results_dict)   
        comm.gather(results_dict_pkl, root=0)

    # STEP 3: Save results 
    if rank == 0:
        # Process gathered results
        combined_results = {}
        for result_pkl in results:
            if not result_pkl:
                continue
            result = pickle.loads(result_pkl) 
            for key, value in result.items():
                if key not in combined_results:
                    combined_results[key] = value
                else:
                    combined_results[key]['SCORE'] += value['SCORE']
                    combined_results[key]['ALLELES_FOUND'] += value['ALLELES_FOUND']
                    combined_results[key]['NUMBER_IMPUTED'] += value['NUMBER_IMPUTED']
                    if config['PGSCalc']['incude_snps_per_sample']:
                        combined_results[key]['SNP'].append(value['SNP'])

        if not combined_results:
            logging.error("Combined results empty... aborting all ranks.")
            comm.Abort()

        logging.info("Combined results DataFrame created.")
        results_df = pd.DataFrame.from_dict(combined_results, orient='index')

        output_filename = generate_output_filename(config)
        
        # Save results for each PGS_ID separately
        for pgs_id in score_metadata.keys():
            pgs_results_df = results_df[results_df['PGS_ID'] == pgs_id]
            if config["PGSCalc"]["summary_statistics"]:
                calculate_and_save_summary_statistic(pgs_results_df, config, output_filename, pgs_id)
            save_results(pgs_results_df, config, output_filename, pgs_id)
            logging.info(f"Results for PGS_ID {pgs_id} saved.")

        logging.info("Completed.")

def sequential_main(config, args):
    """Sequential PGSCalc pipeline."""

    setup_logging(args.job_id)

    # Pull scoring schema
    score_package = get_combined_scores(config)
    num_of_snps = sum(len(list(sub_dict.values())[0]) for sub_dict in score_package.combined_scores.values() if sub_dict.values())

    if num_of_snps == 0:
        logging.error("No SNPs to query. Check if accessions are valid, or attempt PGScore override if files are corrupted.")
    logging.info(f"Beginning query for {sum(len(list(sub_dict.values())) for sub_dict in score_package.combined_scores.values() if sub_dict.values())} number of SNPs")

    # Pull input data
    ref = VCFQuery.VCFQuery(vcf_locations=config['VCFQuery']['input_dir'], override=config['VCFQuery']['override'])

    # Set up output vars
    output_filename = generate_output_filename(config)
    score_metadata = {score.pgs_id: sum(len(positions) for positions in score.snps.values()) for score in score_package.scores}

    # Calculate
    scores_to_query = score_package.combined_scores
    if config['PGSCalc']['calculate_overlap_only']:
        scores_to_query = ref.snps_for_cohort(scores_to_query, samples=config['VCFQuery']['sample_ids'])
    
    combined_results = process_snps(
        scores_to_query, ref, score_metadata, 
        samples=config['VCFQuery']['sample_ids'], 
        remove_imputed=config['VCFQuery']['remove_imputed'], 
        remove_filtered=config['VCFQuery']['remove_filtered'],
        include_snps=config['PGSCalc']['incude_snps_per_sample']
        )

    # Create DataFrame from combined results
    results_df = pd.DataFrame.from_dict(combined_results, orient='index')

    # Save results for each PGS_ID separately
    for pgs_id in score_metadata.keys():
        pgs_results_df = results_df[results_df['PGS_ID'] == pgs_id]
        if config["PGSCalc"]["summary_statistics"]:
            calculate_and_save_summary_statistic(pgs_results_df, config, output_filename, pgs_id)
        save_results(pgs_results_df, config, output_filename, pgs_id)
        logging.info(f"Results for PGS_ID {pgs_id} saved.")

    logging.info("Completed.")
    
def main():
    args = parse_config()
    config = load_config(args.config)
    if config["PGSCalc"]["parallel"]:
        try:
            multiprocessing_main(config,args)
        except ValueError as e:
            logging.error(e)
            #TODO might need to un-comment later
            sequential_main(config,args)
    else: 
        sequential_main(config,args)

if __name__ == "__main__":
    main()
