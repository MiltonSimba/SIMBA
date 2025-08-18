import argparse
from simba.pipeline import run_pipeline

def main():
    parser = argparse.ArgumentParser(description="SIMBA: Selection Inference & Mutation-Based Annotation")
    parser.add_argument("-f", "--fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("-a", "--annotation", required=True, help="Path to domain annotation CSV")
    parser.add_argument("-t", "--threshold", type=float, default=1.5, help="dN/dS threshold for hotspot detection")
    parser.add_argument("-o", "--output", default="simba_output", help="Prefix for output files")

    args = parser.parse_args()
    run_pipeline(args.fasta, args.annotation, args.threshold, args.output)

if __name__ == "__main__":
    main()

