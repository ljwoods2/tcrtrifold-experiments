params.dset_name = null
params.input_pattern = null


process JOIN_TRIAD_FEATURES {
    label "process_local"
    conda "envs/env.yaml"

    publishDir "${params.data_dir}/${params.dset_name}/triad/staged", mode: 'copy'

    input:
    path base_triad_file
    val triad_feature_glob

    output:
    path("*.parquet")

    script:
    """
    #!/usr/bin/env python
    import polars as pl
    from pathlib import Path

    feat_fnames = "${triad_feature_glob.join(',')}".split(',')

    base_triad = pl.read_parquet("${base_triad_file}")
    
    out_triad = base_triad

    for feat_fname in feat_fnames:
        triad_features = pl.read_parquet(feat_fname)
        out_triad = out_triad.join(triad_features, on="job_name")

    out_triad.write_parquet("${base_triad_file.getBaseName(2)}.feat.parquet")
    """
}

process JOIN_PMHC_FEATURES {
    label "process_local"
    conda "envs/env.yaml"

    publishDir "${params.data_dir}/${params.dset_name}/pmhc/staged", mode: 'copy'

    input:
    path base_pmhc_file
    val triad_feature_glob

    output:
    path("*.parquet")

    script:
    """
    #!/usr/bin/env python
    import polars as pl
    from pathlib import Path

    feat_fnames = "${triad_feature_glob.join(',')}".split(',')

    base_pmhc = pl.read_parquet("${base_pmhc_file}")
    
    out_pmhc = base_pmhc

    for feat_fname in feat_fnames:
        pmhc_feat = pl.read_parquet(feat_fname)
        out_pmhc = out_pmhc.join(pmhc_feat, on="job_name")

    out_pmhc.write_parquet("${base_pmhc_file.getBaseName(2)}.feat.parquet")
    """
}

workflow {
    base_triad_file = "${data_dir}/${dset_name}/triad/features/${input_pattern}"

    base_triad_channel = Channel.fromPath(base_triad_file)

    if base_triad_channel.count() != 1 {
        error "Expected exactly one file matching pattern '${base_triad_file}', but found ${base_triad_channel.count()} files."
    }

    triad_feature_glob = Channel.fromPath("${data_dir}/${dset_name}/triad/features/${input_pattern}")

    JOIN_TRIAD_FEATURES(
        base_triad_file,
        triad_feature_glob
    )

    base_pmhc_file = "${data_dir}/${dset_name}/pmhc/features/${input_pattern}"
    base_pmhc_channel = Channel.fromPath(base_pmhc_file)

    if base_pmhc_channel.count() != 1 {
        error "Expected exactly one file matching pattern '${base_pmhc_file}', but found ${base_pmhc_channel.count()} files."
    }

    pmhc_feature_glob = Channel.fromPath("${data_dir}/${dset_name}/pmhc/features/${input_pattern}").collect()

    JOIN_PMHC_FEATURES {
        base_pmhc_file,
        pmhc_feature_glob
    }
}