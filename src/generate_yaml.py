import argparse
from pathlib import Path
import yaml

def generate_yaml(pool, sample_id, tcr='A3'):
    tcr_params = {
        'A3': {
            'wt': "TCAGGGCGCCAGTTCTCTAACTCTCGCTCTGAGATGAATGTGAGCACCTTGGAGCTGGGGGACTCGGCCCTTTATCTTTGCGCCAGCAGCCCGAATATGGCGGATGAACAGTACTTCGGGCCGGGCACCAGGCTCACGGTCACAGGTAAGGCTGGGGGTCTATAGGAGGGGTGCGATGAGGGAGGACTCTGTCCTGGGAAATGTCAAAGAGAACAGAGATCCCAGCTCCCGGAGCCAGACTGAGGGAGACGTCATGT",
            'fwd_adapter': "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
            'rc_adapter': "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
            'expected_length': 269
        }
        # Add more TCRs as needed
    }

    if tcr not in tcr_params:
        raise ValueError(f"TCR {tcr} not recognized. Available: {list(tcr_params.keys())}")

    config = {
        'quality_threshold': 20,
        'min_length': 50,
        'expected_length': 297,
        'seq_threshold': 5,
        'qtrim': "r",
        
        **tcr_params[tcr],

        'repair': True,
        'trim': True,
        'filter': True,

        'pool': pool,
        'sample_id': sample_id,
        'TCR': tcr,
        'input_dir': f'/cluster/project/reddy/marluca/NGS_pipeline/data/raw/{pool}',
        'output_dir': f'/cluster/project/reddy/marluca/NGS_pipeline/data/processed/{pool}',
        
        'bbmap_dir': "/cluster/project/reddy/marluca/NGS_pipeline/bbmap"
    }

    config_dir = Path('/cluster/project/reddy/marluca/NGS_pipeline/config')
    config_dir.mkdir(parents=True, exist_ok=True)
    config_path = config_dir / f"{tcr}_CDR3_library-config_{sample_id}.yaml"

    with open(config_path, 'w') as f:
        yaml.dump(config, f, sort_keys=False)

    print(f"✅ Wrote {config_path} (TCR={tcr})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate YAML config for samples and TCR")
    parser.add_argument('--samples', type=str, nargs='+', required=True, help='List of sample IDs')
    parser.add_argument('--pool', type=str, required=True, help='sample pool')
    parser.add_argument('--TCR', type=str, default='A3', help='TCR type (default: A3)')

    args = parser.parse_args()

    for sample in args.samples:
        generate_yaml(args.pool, sample, args.TCR)
