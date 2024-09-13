import os

def find_vcf_files(input_dir):
    # Initialize a dictionary to store the directory structure
    dir_structure = {}

    # Walk through the directory and its subdirectories
    for root, _, files in os.walk(input_dir):
        # Print the directory being scanned
        print(f"Scanning directory: {root}")

        # Initialize the list of VCF files for the current directory
        vcf_files = []

        # Filter for files ending in .gvcf, .gvcf.gz, .vcf, .vcf.gz
        for f in files:
            if f.endswith(('.gvcf', '.gvcf.gz', '.vcf', '.vcf.gz')):
                print(f"Found file: {f}")
                vcf_files.append(f)
        
        if vcf_files:
            # Add the directory and its VCF files to the structure
            relative_dir = os.path.relpath(root, input_dir)
            dir_structure[relative_dir] = vcf_files

    return dir_structure

def print_directory_structure(dir_structure):
    print(f"{input_dir}:")
    for subdir, files in dir_structure.items():
        print(f"|-----{subdir}: {', '.join(files)}")

if __name__ == "__main__":
    # Replace with your input directory
    input_dir = "/N/project/biobank/phi_ingest_biobank_regeneron_KEEPTRXyunlong"

    # Find the VCF files and print the directory structure
    dir_structure = find_vcf_files(input_dir)
    print_directory_structure(dir_structure)
