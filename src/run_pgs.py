import argparse
import json
import hashlib
import os
from PGScore import PGScore
import subprocess

# TODO:
#   -[X] Add ability for multithreading for large PGS objects (most should use it automatically!) but we want to scale accordingly 
#   -[X] When using trait ID, will need to make API call to identify unique scores 
#   -[] If usinng trait ID, want some organization to store results together ... 

def parse_config():
    """Parse the configuration file argument."""
    parser = argparse.ArgumentParser(description="Download PGS scoring files and query sample data.")
    parser.add_argument('-c', '--config', required=True, help="Path to the configuration file.")
    parser.add_argument('-s', '--slurm', required=True, help='Path to the template Slurm script.')
    parser.add_argument('-cdir', '--config_dir', required=True, help="Directory for individual score configs to be saved.")
    parser.add_argument('-sdir', '--script_dir', required=True, help="Directory for individual score scripts to be saved.")
    return parser.parse_args()

def load_config(config_path):
    """Load the configuration file."""
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config

def get_score_ids(config):
    """Pull score IDs either using config or PGScore."""
    score_package = PGScore.ScorePackage(
        pgs_ids=config['PGScore'].get('pgs_ids', None), 
        trait_id=config['PGScore'].get('trait_id', None),
        trait_include_children=config['PGScore'].get('trait_include_children', 1), 
        genome_build=config['PGScore']['genome_build'],
        override=config['PGScore']['override']
    )
    ids_nodes = []
    for score in score_package.scores:
        # print(score.pgs_id)
        nodes = len(score.snps) // 500000
        ceiling = 20
        ids_nodes.append([score.pgs_id, min(nodes, ceiling)])
    return ids_nodes


def generate_config(original_config, score, c_destination):
    new_config = original_config
    new_config["PGScore"]["pgs_ids"] = [score]

    new_config["PGScore"].pop('trait_id',0)
    new_config["PGScore"].pop('trait_include_children',0)

    unique_string = hashlib.sha256(json.dumps(new_config, sort_keys=True).encode('utf-8')).hexdigest()[:8]
    config_fp = os.path.join(c_destination, f'{unique_string}.json')
    with open(config_fp, 'w') as f:
        json.dump(new_config, f, indent=4)
    return config_fp

def generate_script(slurm_template, config_fp, nodes, s_destination):
    with open(slurm_template, 'r') as file:
        lines = file.readlines()

    unique_string = os.path.splitext(os.path.basename(config_fp))[0]
    script_fp = os.path.join(s_destination, f'{unique_string}.script')

    with open(script_fp, 'w') as file:
        for line in lines:
            if '#SBATCH --nodes=' in line:
                current_nodes = int(line.split('=')[1].strip())
                new_nodes = max(current_nodes, nodes)
                # line = f'#SBATCH --nodes={new_nodes}\n'
                line = f'#SBATCH --nodes=1\n'

            if 'srun python' in line and '--config' in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if '--config' in part:
                        parts[i + 1] = config_fp
                        break
                line = ' '.join(parts) + '\n'

            file.write(line)

    return script_fp 

def execute_script(script_fp):
    # print(f"Executing: {script_fp}")
    subprocess.run([
        'sbatch', script_fp
    ], check=True)

def main():
    args = parse_config()
    config = load_config(args.config)
    slurm_template = args.slurm
    for score, nodes in get_score_ids(config):
        print(f"Processing score: {score} with {nodes} nodes")
        
        # Generate config file and print the path
        score_config_fp = generate_config(config, score, args.config_dir)
        print(f"Generated config file at: {score_config_fp}")
        
        # Generate script file and print the path
        script_fp = generate_script(slurm_template, score_config_fp, nodes, args.script_dir)
        print(f"Generated script file at: {script_fp}")
        
        # Execute the script and print confirmation
        execute_script(script_fp)
        print(f"Executed script for score {score}\n")

if __name__ == "__main__":
    main()
