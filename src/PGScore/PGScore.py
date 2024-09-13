import subprocess
import os
import pandas as pd
import logging
import pickle
import hashlib
import requests 

# Get the logger instance
logger = logging.getLogger(__name__)

class ScorePackage:
    def __init__(self, pgs_ids=None, trait_id=None, trait_include_children=False, genome_build='GRCh38', scoring_file_directory=None, override=True):
        logging.info("Initializing score package.")
        logging.info(f"trait: {trait_id}")
        self.output_dir = "./src/PGScore/scoring_files/"
        self.scores = []
        self.combined_scores = {}

        if scoring_file_directory:
            self._load_local_scores(scoring_file_directory)
        else:
            if pgs_ids:
                subset_hash = hashlib.md5(','.join(sorted(pgs_ids)).encode()).hexdigest()[:8]
            elif trait_id:
                subset_hash = hashlib.md5(trait_id.encode()).hexdigest()[:8]

            if subset_hash:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                save_path = os.path.join(script_dir, "../PGScore/score_package_objects", f"{subset_hash}_{genome_build}_ref.pkl")
                save_path = os.path.abspath(save_path)
            else:
                logging.warning("Unable to generate unique ID for subset. ScorePackage not cacheable.")

            if not override and os.path.exists(save_path):
                self.load(save_path)
                logging.info("Using previously saved ScorePackage object.")
                
            else:
                self.pgs_ids = pgs_ids # Don't know this if using PGS ID
                self.trait_id = trait_id
                self.trait_include_children = trait_include_children
                self.genome_build = genome_build
                self._download_scores(override)
                self._combine_scores()
                if save_path:
                    self.save(save_path)
            logging.info(f"Scoring object built for ids: {self.pgs_ids}")

    def save(self, save_path):
        with open(save_path, 'wb') as f:
            pickle.dump(self, f)
        logging.info(f"VCFQuery object saved to {save_path}.")

    def load(self, save_path):
        with open(save_path, 'rb') as f:
            obj = pickle.load(f)
        self.__dict__.update(obj.__dict__)
        logging.info(f"VCFQuery object loaded from {save_path}.")

    def _download_scores(self, override):
        logging.info("Starting the download of PGS scores.")

        if self.trait_id:
            logging.info(f"Pulling associated traits with {self.trait_id}")
            base_url = "https://www.pgscatalog.org/rest/trait/"
            url = f"{base_url}{self.trait_id}?include_children={self.trait_include_children}"
            response = requests.get(url)
            
            if response.status_code == 200:
                data = response.json()
                pgs_ids = data.get("associated_pgs_ids", [])
                logging.info(f"Received PGS IDs: {pgs_ids}")
                # TODO if we want to change so you can do a trait and input pgs ids, change here!
                self.pgs_ids = list(pgs_ids)
            
            command = f"python -m pgscatalog_utils.download.download_scorefile -t {self.trait_id} -b {self.genome_build} -o {self.output_dir}"
        else:
            ids_str = ' '.join(self.pgs_ids)
            command = f"python -m pgscatalog_utils.download.download_scorefile -i {ids_str} -o {self.output_dir} -b {self.genome_build}"
        # Run the command
        try:
            subprocess.run(command, shell=True, check=True)
            logging.info("Download completed successfully.")
        except subprocess.CalledProcessError as e:
            logging.error(f"An error occurred while downloading scores: {e}")
            return

        # Create PGScore objects
        for pgs_id in self.pgs_ids:
            # TODO double check this - harmonized is getting fetched
            scoring_file_path = os.path.join(self.output_dir, f"{pgs_id}_hmPOS_{self.genome_build}.txt.gz")
            if os.path.exists(scoring_file_path):
                logging.info(f"Creating PGScore object for {pgs_id}.")
                try:
                    self.scores.append(PGScore(pgs_id, scoring_file_path, self.genome_build, override))
                except Exception as e:
                    logging.error(f"Error creating PGScore object for {pgs_id}: {e}")
            else:
                logging.warning(f"Warning: {scoring_file_path} does not exist.")

    def _load_local_scores(self, directory):
        logging.info("Loading local scorefiles.")
        for filename in os.listdir(directory):
            if filename.endswith('.txt') or filename.endswith('.txt.gz'):
                pgs_id = filename.split('_')[0]
                scoring_file_path = os.path.join(directory, filename)
                if os.path.exists(scoring_file_path):
                    logging.info(f"Creating PGScore object for {pgs_id}.")
                    self.scores.append(PGScore(pgs_id, scoring_file_path))
                else:
                    logging.warning(f"Warning: {scoring_file_path} does not exist.")

    def _combine_scores(self):
        for score in self.scores:
            snp = score.snps
            # print(f"adding {snp}")
            for chrom, pos_dict in snp.items():
                if chrom not in self.combined_scores.keys():
                    self.combined_scores[chrom] = {}
                for pos, scoring_info in pos_dict.items():
                    if pos not in self.combined_scores[chrom].keys():
                        self.combined_scores[chrom][pos] = []
                    self.combined_scores[chrom][pos].extend(snp[chrom][pos])

class PGScore:
    def __init__(self, pgs_id, scoring_file_path, build="GRCh38", override=True):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        save_path = os.path.join(script_dir, "../PGScore/score_objects", f"{pgs_id}_{build}_ref.pkl")
        save_path = os.path.abspath(save_path)

        if not override and os.path.exists(save_path):
            self.load(save_path)
            logging.info(f"Using previously saved PGScore object for {pgs_id} in {build}.")
        else:
            self.pgs_id = pgs_id
            self.scoring_file_path = scoring_file_path
            self.data = self._parse_scoring()
            self.snps = self._extract_snps()
            self.save(save_path)

    def _parse_scoring(self):
        logging.info(f"Parsing scoring file for {self.pgs_id}.")
        try:
            if self.scoring_file_path.endswith('.gz'):
                data = pd.read_csv(self.scoring_file_path, compression='gzip', sep='\t', comment='#')
            else:
                data = pd.read_csv(self.scoring_file_path, sep='\t', comment='#')
        except Exception as e:
            logging.error(f"Error reading scoring file {self.scoring_file_path}: {e}")
            raise
        return data

    def _extract_snps(self):
        logging.info(f"Extracting SNPs from scoring file for {self.pgs_id}.")
        snps = {}
        haplotype_skipped = False
        # TODO need more checking for complex data types
        for _, row in self.data.iterrows():
            if 'is_haplotype' in row and row['is_haplotype']:
                haplotype_skipped = True
                continue

            try:
                try:
                    snp_chrom = str(int(row['hm_chr']))
                except ValueError:
                    if row['hm_chr'] in ['X', 'Y']:
                        snp_chrom = row['hm_chr']
                    else:
                        raise ValueError(f"Unexpected chromosome value: {row['hm_chr']}")
                snp_pos = int(row['hm_pos'])

                if snp_chrom not in snps.keys():
                    snps[snp_chrom] = {}
                if snp_pos not in snps[snp_chrom].keys():
                    snps[snp_chrom][snp_pos] = []

                snps[snp_chrom][snp_pos].extend([{'score': self.pgs_id, 'effect_allele': row['effect_allele'],'other_allele': tuple(row['other_allele'].split(',')) if 'other_allele' in row else None,'rsID': row['hm_rsID'],'effect_weight': row['effect_weight']}])

            except Exception as e: 
                logging.warning(f"For PGS {self.pgs_id}, due to error {e} unable to add row {row}.")

        if haplotype_skipped:
            warning_message = f"Warning: Haplotype SNPs are skipped in PGS ID {self.pgs_id}. Using this scoring file may result in an inaccurate score."
            logging.warning(warning_message)

        return snps

    def save(self, save_path):
        with open(save_path, 'wb') as f:
            pickle.dump(self, f)
        logging.info(f"PGScore object saved to {save_path}.")

    def load(self, load_path):
        with open(load_path, 'rb') as f:
            loaded_obj = pickle.load(f)
        self.__dict__.update(loaded_obj.__dict__)
        logging.info(f"PGScore object loaded from {load_path}.")

# class PGScore:
#     def __init__(self, pgs_id, scoring_file_path):
#         self.pgs_id = pgs_id
#         self.scoring_file_path = scoring_file_path
#         self.data = self._parse_scoring()
#         self.snps = self._extract_snps()

#     def _parse_scoring(self):
#         logging.info(f"Parsing scoring file for {self.pgs_id}.")
#         try:
#             if self.scoring_file_path.endswith('.gz'):
#                 data = pd.read_csv(self.scoring_file_path, compression='gzip', sep='\t', comment='#')
#             else:
#                 data = pd.read_csv(self.scoring_file_path, sep='\t', comment='#')
#         except Exception as e:
#             logging.error(f"Error reading scoring file {self.scoring_file_path}: {e}")
#             raise
#         return data
        
#     def _extract_snps(self):
#         logging.info(f"Extracting SNPs from scoring file for {self.pgs_id}.")
#         snps = []
#         haplotype_skipped = False
#         for _, row in self.data.iterrows():
#             if 'is_haplotype' in row and row['is_haplotype']:
#                 haplotype_skipped = True
#                 continue
#             snp = {
#                 'chrom': str(row['hm_chr']),
#                 'pos': row['hm_pos'],
#                 'effect_allele': row['effect_allele'],
#                 'other_allele': tuple(row['other_allele'].split(',')),
#                 'rsID': row['hm_rsID'],
#                 'effect_weight': row['effect_weight']
#             }
#             snps.append(snp)
#         if haplotype_skipped:
#             warning_message = f"Warning: Haplotype SNPs are skipped in PGS ID {self.pgs_id}. Using this scoring file may result in an inaccurate score."
#             # print(warning_message)
#             logging.warning(warning_message)
#         return snps

# if __name__ == "__main__":
#     # Example usage
#     pgs_ids = ['PGS000001']
#     output_dir = './scoring_files/'
#     downloader = PGSDownloader(pgs_ids, output_dir)
#     downloader.download_scores()

#     for score in downloader.scores:
#         logging.info(f"PGS ID: {score.pgs_id}")
#         logging.info(f"First 5 SNPs: {score.snps[0:5]}")
#         print(len(score.snps))
