#!/bin/bash
#SBATCH --job-name=blastpdb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwoods@tgen.org
#SBATCH --ntasks=1
#SBATCH --mem=3G
#SBATCH --time=72:00:00
#SBATCH --output=run_alpha.%j.log

. ~/miniconda3/etc/profile.d/conda.sh
conda activate blast

blastp -query fasta_queries/all_triads.fasta -db /tgen_labs/altin/alphafold3/blast/pdbaa -out pdb_blast_results/blast_result.csv -outfmt 10
