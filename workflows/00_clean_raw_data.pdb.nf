"""
Hardcoded pipelines for cleaning raw data files from each dataset, getting them into a uniform format.
"""

process CLEAN_PDB {
  label "process_local"
  conda "envs/env.yaml"

  publishDir(
    path: {"${params.data_dir}/pdb/triad/staged"},
    pattern: "pdb.triad*",
    mode: 'copy'
  )
  publishDir(
    path: {"${params.data_dir}/pdb/pmhc/staged"},
    pattern: "pdb.pmhc*",
    mode: 'copy'
  )


  input:
  path pdb_rep
  path pdb_stcr

  output:
  path("*.parquet")

  script:
  """
  clean_pdb.py \
    --raw_csv_path ${pdb_rep} \\
    --raw_stcr_path ${pdb_stcr} \\
    --imgt_hla_path ${params.imgt_hla_path} \\
    -ot pdb.triad.cleaned.parquet \\
    -op pdb.pmhc.cleaned.parquet 
  """
}


workflow {

  CLEAN_PDB(Channel.fromPath("data/pdb/raw/table_S1_structure_benchmark_complexes.csv"),
    Channel.fromPath("data/pdb/raw/db_summary.dat"))

}