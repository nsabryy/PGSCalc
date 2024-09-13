import os 
import pysam
import logging
import re
import pickle

# Get the logger instance
logger = logging.getLogger(__name__)

class VCFSource(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.filename = os.path.basename(filepath)  # store filename 
        self.samples = set()  # store samples included
        self.chromosomes = set()  # store chromosomes present
        self.build = None  # store genome build
        self._parse_vcf_header()
        self.index_filename = None

    def _parse_vcf_header(self):
        try:
            if os.path.exists(self.filepath + '.csi'):
                vcf_in = pysam.VariantFile(self.filepath, 'r', index_filename=self.filepath + '.csi')
                self.index_filename = self.filepath + '.csi'
            elif os.path.exists(self.filepath + '.tbi'):
                vcf_in = pysam.VariantFile(self.filepath, 'r', index_filename=self.filepath + '.tbi')
                self.index_filename = self.filepath + '.tbi'
            else:
                try: 
                    vcf_in = pysam.VariantFile(self.filepath, 'r')
                except:
                    raise FileNotFoundError("No index file found for " + self.filepath)

            self.samples = set(vcf_in.header.samples)
            self.chromosomes = set(contig.lstrip("chr") for contig in vcf_in.header.contigs)
            if not self.chromosomes:
                self.chromosomes = set(str(i) for i in range(1, 23))

            self.build = self._extract_genome_build(vcf_in.header)

        except (OSError, IOError, pysam.libcbcf.BCFError) as e:
            logging.error(f"Error reading VCF file {self.filepath}: {e}")
            raise  # re-raise the exception to handle it in the VCFQuery class

    def _extract_genome_build(self, header):
        for line in header.records:
            if line.key == 'reference':
                reference = line.value.lower()
                if 'grch37' in reference or 'hg19' in reference:
                    return 'GRCh37'
        return 'GRCh38'

    def check_sampleid(self, check_id):
        return (check_id in self.samples)
    
    def check_chromosomes(self, check_chromosome):
        return (check_chromosome in self.chromosomes)

class VCFQuery:
    def __init__(self, vcf_locations=[], override=False):
        self._vcfs = []  # list of possible VCF genotype sources
        self._samples = {}  # map sample IDs to VCF files
        self._chromosomes = {} # map chromosomes to VCF files
        self.inaccessible_files = 0

        # Get full path for VCFQuery pickles
        script_dir = os.path.dirname(os.path.abspath(__file__)) 

        for vcf_location in vcf_locations:
            # Normalize the directory path for saving
            input_dir_normalized = re.sub(r'[\\/:*?"<>|]', '_', vcf_location.strip('/').replace(' ', '_'))
            save_path = os.path.join(script_dir, "../VCFQuery/query_objects", f"{input_dir_normalized}_ref.pkl")
            save_path = os.path.abspath(save_path)

            # Load a saved VCFDirectory object if available and not overridden
            if not override and os.path.exists(save_path):
                self.load(save_path)
                logging.info(f"Using previously saved VCFDirectory object for {vcf_location}.")
            else:
                # Add sources from the directory
                self.add_source(vcf_location)
                
                # Save the object after processing the directory
                if save_path:
                    self.save(save_path)

    def get_save_path(self, vcf_location):
        script_dir = os.path.dirname(os.path.abspath(__file__)) 
        input_dir_normalized = re.sub(r'[\\/:*?"<>|]', '_', vcf_location.strip('/').replace(' ', '_'))
        save_path = os.path.join(script_dir, "../VCFQuery/query_objects", f"{input_dir_normalized}_ref.pkl")
        save_path = os.path.abspath(save_path)
        return save_path

    def save(self, save_path):
        with open(save_path, 'wb') as f:
            pickle.dump(self, f)
        logging.info(f"VCFQuery object saved to {save_path}.")

    def load(self, save_path):
        with open(save_path, 'rb') as f:
            obj = pickle.load(f)
        
        self._vcfs.extend(obj._vcfs)
        self._samples.update(obj._samples)
        self._chromosomes.update(obj._chromosomes)
        self.inaccessible_files += obj.inaccessible_files
        
        logging.info(f"VCFQuery object loaded and merged from {save_path}.")

    def add_source(self, src):
        if os.path.isdir(src):
            # Check if the directory is readable
            if os.access(src, os.R_OK):
                # find vcfs in path
                print(f"Scanning directory: {src}")
                for root, dirs, files in os.walk(src):
                    for file in files:
                        if file.endswith('.vcf') or file.endswith('.vcf.gz') or file.endswith('.gvcf') or file.endswith('.gvcf.gz'):
                            file_path = os.path.join(root, file)
                            if os.access(file_path, os.R_OK):
                                self.add_source(file_path)
                            else:
                                self.inaccessible_files += 1
            else:
                self.inaccessible_files += 1
                logging.error(f"Permission denied for directory {src}")
        else:
            # Check if the file is readable and TBI/CSI index exists
            #  and ((os.path.exists(src + '.tbi') or os.path.exists(src + '.csi')))
            if os.access(src, os.R_OK):
                try:
                    vcf_source = VCFSource(src)
                    self._vcfs.append(vcf_source)

                    # populate samples
                    for sample in vcf_source.samples:
                        if sample not in self._samples:
                            self._samples[sample] = []
                        self._samples[sample].append(vcf_source.filename)

                    for chromosome in vcf_source.chromosomes:
                        if chromosome not in self._chromosomes:
                            self._chromosomes[chromosome] = []
                        self._chromosomes[chromosome].append(vcf_source.filename)

                except Exception:
                    self.inaccessible_files += 1
                    logging.error(f"Error adding VCF source {src}")
            else:
                self.inaccessible_files += 1
                logging.error(f"Permission denied or index missing for file {src}")
    
    def query_sources(self, chrom=[], samples=[]):
        chrom_set = set()
        if chrom:
            for c in chrom:
                if c in self._chromosomes.keys():
                    chrom_set.update(self._chromosomes[c])

        sample_set = set()
        if samples:
            for s in samples:
                if s in self._samples.keys():
                    sample_set.update(self._samples[s])

        if not sample_set and not chrom_set: 
            result_filenames = set(self._vcfs)
        elif not sample_set or not chrom_set:
            result_filenames = chrom_set.union(sample_set)
        else:
            result_filenames = chrom_set.intersection(sample_set)
        
        result_sources = [vcf for vcf in self._vcfs if vcf.filename in result_filenames]
        return result_sources
    
    def _format_call(self, call, ref, alts):
        alts = alts if alts is not None else []
        allele_list = [ref] + list(alts)
        formatted_call = tuple(
            allele_list[c] if c is not None and c != '.' else c for c in call
        )
        return formatted_call

    def _check_imputed(self, record):
        return 'IMPUTED' in record.info

    def _check_filter(self, record):
        return record.filter.keys() not in [[], ['PASS']]

    def _allele_verification(self, ref, alt, eff, oth):
        if oth == set(['.']):
            if not (
                (eff == ref) or 
                (eff == alt)
            ):
                return False
        else:
            if not (
                (ref == eff and alt == oth) or 
                (ref == oth and alt == eff)
            ):
                return False
        return True

    def fetch_alleles(self, chrom, pos):
        try:
            for vcf_source in self.query_sources(chrom=[chrom]):
                vcf_in = pysam.VariantFile(vcf_source.filepath, 'r')
                for record in vcf_in.fetch(chrom, pos - 1, pos):
                    return record.ref, record.alts
        except Exception as e:
            logging.warning(f"Failed to fetch alleles for {chrom}:{pos} due to {e}")
        return None, None

    def query_snps(self, snps, samples=None, remove_imputed=False, remove_filtered=False, only_imputed=False, overlap=True):
        # Verify user filters.
        if remove_imputed == True and only_imputed==True:
            logging.error(f"remove_imputed and only_imputed both True meaning no snps will qualify. Returning None")
            return
            
        # Pull relevant VCF files to query.
        vcf_sources_by_chrom = {}
        for chrom in snps.keys():
            vcf_sources_by_chrom[chrom] = self.query_sources(chrom=[chrom], samples=samples)
            
            if not vcf_sources_by_chrom[chrom]:
                logging.warning(f"No matching VCF sources found for chromosome {chrom}")

        # Format samples if present. 
        if samples: 
            samples = set(samples)

        # Begin SNP query. 
        for chrom, pos_list in snps.items():
            # Check chromosome data exists. 
            if chrom not in vcf_sources_by_chrom or not vcf_sources_by_chrom[chrom]:
                continue
            chrom_query = chrom
            # Iterate through given positions. 
            for pos in pos_list:    
                snp_in_vcf = False
                for vcf_source in vcf_sources_by_chrom[chrom]:
                    try:
                        if vcf_source.index_filename:
                            vcf_in = pysam.VariantFile(vcf_source.filepath, 'r', index_filename=vcf_source.index_filename)
                        else: 
                            vcf_in = pysam.VariantFile(vcf_source.filepath, 'r')
                        try:
                            records = vcf_in.fetch(chrom, pos - 1, pos)

                        except:
                            chrom_with_prefix = f"chr{chrom}"
                            records = vcf_in.fetch(chrom_with_prefix, pos - 1, pos)

                        for record in records:
                            if record.pos == pos:
                                # When SNP found, perform requested filtering upon records. 
                                snp_in_vcf = True

                                if record.alts is not None and any(len(a) > 1 for a in record.alts):
                                    logging.warning(f"In file {vcf_source.filename}: MNP encountered at {chrom}:{pos}. REF={record.ref} ALT={record.alts}. Logic for MNP is not set up yet.")
                                    continue

                                if remove_imputed and self._check_imputed(record):
                                    logging.info(f"In file {vcf_source.filename}: Skipped imputed SNP: {chrom}:{pos} REF={record.ref} ALT={record.alts}")
                                    continue

                                if only_imputed and not self._check_imputed(record):
                                    continue

                                if remove_filtered and self._check_filter(record):
                                    logging.info(f"In file {vcf_source.filename}: Skipped filtered SNP: {chrom}:{pos} REF={record.ref} ALT={record.alts} FILTER={record.filter.keys()}")
                                    continue

                                # Yield SNP record to PGSCalc. 
                                yield vcf_source.filepath, record

                    except (OSError, IOError) as e:
                        logging.error(f"Error reading VCF file {vcf_source.filepath}: {e}")
                
                if not snp_in_vcf:
                    logging.warning(f"No matching VCF records for {chrom}:{pos}.")

    def snps_for_cohort(self, combined_scores, samples=None):
        snps = {chrom: [int(key) for key in pos_dict.keys()] for chrom, pos_dict in combined_scores.items()}
        if samples:
            all_samples = set(samples)
        else:
            all_samples = set(self._samples.keys())
        vcf_sources_by_chrom = {}
        for chrom in snps.keys():
            vcf_sources_by_chrom[chrom] = self.query_sources(chrom=[chrom])

        for chrom, pos_list in snps.items():
            if chrom not in vcf_sources_by_chrom or not vcf_sources_by_chrom[chrom]:
                continue
            chrom_query = chrom
            for pos in pos_list:    
                temp_list = all_samples.copy()
                for vcf_source in vcf_sources_by_chrom[chrom]:
                    try:
                        if vcf_source.index_filename:
                            vcf_in = pysam.VariantFile(vcf_source.filepath, 'r', index_filename=vcf_source.index_filename)
                        else: 
                            vcf_in = pysam.VariantFile(vcf_source.filepath, 'r')                        
                        try:
                            records = vcf_in.fetch(chrom, pos - 1, pos)

                        except:
                            chrom_with_prefix = f"chr{chrom}"
                            records = vcf_in.fetch(chrom_with_prefix, pos - 1, pos)

                        for record in records:
                            if record.pos == pos:
                                temp_list -= vcf_source.samples
                    except (OSError, IOError) as e:
                        logging.error(f"Error reading VCF file {vcf_source.filepath}: {e}")
                if temp_list:
                    combined_scores[chrom].pop(pos)
                    logging.info(f"Removed {chrom}:{pos} from scoring object due to missing sample data.")
            if not combined_scores[chrom]:
                combined_scores.pop(chrom)
        return combined_scores
