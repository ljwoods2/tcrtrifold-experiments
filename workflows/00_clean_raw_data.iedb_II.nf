"""
Hardcoded pipelines for cleaning raw data files from each dataset, getting them into a uniform format.
"""

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
  clean_iedb_II.py \\
    --raw_csv_path ${iedb_II} \\
    -ot iedb_II.triad.cleaned.parquet \\
    -op iedb_II.pmhc.cleaned.parquet 
  """
}


workflow {

  CLEAN_IEDB_II(Channel.fromPath("data/iedb_II/raw/immrep_IEDB.csv"))

}