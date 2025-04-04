
## Extracting raw cognate triads from IEDB and VDJDB

The raw extracted data from IEDB and VDJDB is already available in this repo [here](data/iedb-vdjdb/raw/)

Alternatively, if you're interested in re-running our data extraction, clone our fork of the IEDB_IMMREP data repo and run the data extraction script.
```bash
git clone https://github.com/ljwoods2/IEDB_IMMREP.git
cd IEDB_IMMREP
git checkout new-categories
python setup.py install
chmod +x run.sh
./run.sh
```

## Processing cognate triads from IEDB and VDJDB

Processed triads are available by category (species and MHC class) in [data/iedb-vdjdb/iedb](data/iedb-vdjdb/iedb) and [data/iedb-vdjdb/vdjdb](data/iedb-vdjdb/vdjdb/) in the format described in [tcr_format_parsers](https://github.com/ljwoods2/tcr_format_parsers).

If you're interested in re-running our formatting code, first create a conda environment containing the necessary dependencies:

```bash
conda create -n af3-analyzer --file envs/af3-analyzer.yaml
```

Then, run each cell sequentially in [data/iedb-vdjdb/reformat.ipynb](data/iedb-vdjdb/reformat.ipynb) using the `af3-analyzer` environment as the kernel.


## Running AF3 inference for triads from IEDB and VDJDB

Our lab used the nextflow pipelines in the [af3-nf](https://github.com/ljwoods2/af3-nf) repo for running inference on triads. While these pipelines were designed to run on TGen's Gemini supercomputer, they can be easily adapted to run in other environments. Please contact the authors for details.

See [data/iedb-vdjdb/iedb/human_I/run_af3_triad.sh](data/iedb-vdjdb/iedb/human_I/run_af3_triad.sh) for an example slurm script that runs the pipelines.

## Identifying IEDB and VDJDB triad overlap with PDB

Blast+ for pdb alignment

```bash
conda create -n blast --file envs/blast.yaml
```

```bash
cd /path/to/pdbaa/dir
update_blastdb.pl --decompress pdbaa
```

```bash
cd data/iedb-vdjdb
blastp -query fasta_queries/all_triads.fasta -db /path/to/pdbaa/dir -out pdb_blast_results/blast_result.csv -outfmt 10
```

## Processing cognate triads from PDB

Processed PDB triads are available in [data/pdb/pdb_triads.csv](data/pdb/pdb_triads.csv).

Raw PDB summary files in [data/pdb/raw](data/pdb/raw) come from [STCRDab](https://opig.stats.ox.ac.uk/webapps/stcrdab-stcrpred).

If you're interested in re-running our formatting code, first clone the IMGTHLA repo (this is used to identify likely MHC alleles for each sequence):
```bash
git clone https://github.com/ANHIG/IMGTHLA
```

Then, run each cell sequentially in [data/pdb/reformat.ipynb](data/pdb/reformat.ipynb) using the `af3-analyzer` environment as the kernel, making sure to modify the variable `IMGT_HLA_PATH` with your own path to the cloned IMGTHLA repo.


