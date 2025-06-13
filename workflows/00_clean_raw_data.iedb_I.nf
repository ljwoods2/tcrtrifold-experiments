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
  clean_iedb_I.py \\
    --raw_csv_path ${iedb_I} \\
    -ot iedb_I.triad.cleaned.parquet \\
    -op iedb_I.pmhc.cleaned.parquet 
  """
}


workflow {
  CLEAN_IEDB_I(Channel.fromPath("data/iedb_I/raw/immrep_IEDB.csv"))
}