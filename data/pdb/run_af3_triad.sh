#!/bin/bash
#SBATCH --job-name=alphafold
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=3G
#SBATCH --time=72:00:00
#SBATCH --output=run_alpha.%j.log

module load Java/11.0.2

nextflow run \
    -w /scratch/lwoods/work \
    -c /home/lwoods/workspace/af3-nf/nextflow.config \
    /home/lwoods/workspace/af3-nf/af3_triad_msa.nf \
        --input_csv 'pdb_triads.csv' \
        --out_dir '.' \
        --no_peptide \
        -resume \
        --msa_db 'https://pub-vscratch.vast.rc.tgen.org' && \
nextflow run \
    -w /scratch/lwoods/work \
    -c /home/lwoods/workspace/af3-nf/nextflow.config \
    /home/lwoods/workspace/af3-nf/af3_triad_inference.nf \
        --input_csv 'pdb_triads.csv' \
        --out_dir '.' \
        --msa_db 'https://pub-vscratch.vast.rc.tgen.org' \
        --no_peptide \
        --seeds 1,2,3,4,5 \
        --compress \
        -resume
