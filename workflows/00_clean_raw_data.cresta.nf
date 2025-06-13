"""
Hardcoded pipelines for cleaning raw data files from each dataset, getting them into a uniform format.
"""

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
  clean_cresta.py \\
    --raw_csv_path ${cresta} \\
    -ot cresta.triad.cleaned.parquet \\
    -op cresta.pmhc.cleaned.parquet 
  """
}


workflow {

  CLEAN_CRESTA(Channel.fromPath("data/cresta/raw/cresta.csv"))
}