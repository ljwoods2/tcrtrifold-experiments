"""
Hardcoded pipelines for cleaning raw data files from each dataset, getting them into a uniform format.
"""

process CLEAN_IEDB_I {
  label "process_local"
  conda "envs/env.yaml"

  publishDir(
      path: {"${params.data_dir}/iedb_I/triad/staged"},
      pattern: "iedb_I.triad*",
      mode: 'copy'
  )
  publishDir(
      path: {"${params.data_dir}/iedb_I/pmhc/staged"},
      pattern: "iedb_I.pmhc*",
      mode: 'copy'
  )

  input:
  path iedb_I

  output:
  path("*.parquet")

  script:
  """
  clean_iedb_I.py \
    --raw_csv_path ${iedb_I} \
    -ot iedb_I.triad.cleaned.parquet \
    -op iedb_I.pmhc.cleaned.parquet 
  """
}

process CLEAN_IEDB_II {
  label "process_local"
  conda "envs/env.yaml"

  publishDir(
    path: {"${params.data_dir}/iedb_II/triad/staged"},
    pattern: "iedb_II.triad*",
    mode: 'copy'
  )
  publishDir(
    path: {"${params.data_dir}/iedb_II/pmhc/staged"},
    pattern: "iedb_II.pmhc*",
    mode: 'copy'
  )


  input:
  path iedb_II

  output:
  path("*.parquet")

  script:
  """
  clean_iedb_II.py \
    --raw_csv_path ${iedb_II} \
    -ot iedb_II.triad.cleaned.parquet \
    -op iedb_II.pmhc.cleaned.parquet 
  """
}

process CLEAN_CRESTA {
  label "process_local"
  conda "envs/env.yaml"

  publishDir(
    path: {"${params.data_dir}/cresta/triad/staged"},
    pattern: "cresta.triad*",
    mode: 'copy'
  )
  publishDir(
    path: {"${params.data_dir}/cresta/pmhc/staged"},
    pattern: "cresta.pmhc*",
    mode: 'copy'
  )


  input:
  path cresta

  output:
  path("*.parquet")

  script:
  """
  clean_cresta.py \
    --raw_csv_path ${cresta} \
    -ot cresta.triad.cleaned.parquet \
    -op cresta.pmhc.cleaned.parquet 
  """
}

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
    --raw_csv_path ${pdb_rep} \ 
    --raw_stcr_path ${pdb_stcr} \
    -ot pdb.triad.cleaned.parquet \
    -op pdb.pmhc.cleaned.parquet 
  """
}


workflow {

  CLEAN_IEDB_I(Channel.fromPath("data/iedb_I/raw/immrep_IEDB.csv"))
  CLEAN_IEDB_II(Channel.fromPath("data/iedb_II/raw/immrep_IEDB.csv"))
  CLEAN_CRESTA(Channel.fromPath("data/cresta/raw/cresta.csv"))
  CLEAN_PDB(Channel.fromPath("data/pdb/raw/table_S1_structure_benchmark_complexes.csv"),
    Channel.fromPath("data/pdb/raw/db_summary.dat"))  
  // CLEAN_IMMREP(Channel.fromPath("data/immrep/raw/immrep.csv"))
}