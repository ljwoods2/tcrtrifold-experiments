#!/usr/bin/env python
from tcrtrifold.utils import generate_job_name, update_df_from_k_v
from tcrtrifold.pdb import *
from tcr_format_parsers.common.TriadUtils import (
    FORMAT_COLS,
    FORMAT_ANTIGEN_COLS,
    TCRDIST_COLS,
)
from tcr_format_parsers.common.TCRUtils import (
    extract_tcrdist_cols,
)
from mdaf3.FeatureExtraction import serial_apply
import polars as pl
import argparse
from datetime import datetime, timezone


# these serve as a log of changes we make to the PDB data in order to get it into our format
pre_fasta_corrections = {
    "3tf7": {"mhc_chain1": "E", "Achain": None, "Bchain": None},
    # antigen chain located on MHC 2 / TCR 2 chain
    # when we later look for the peptide we can find it by aligning to this chain
    "6bga": {"antigen_chain": "B"},
    "3pl6": {"antigen_chain": "D"},
    "3o6f": {"antigen_chain": "B"},
    "6dfw": {"antigen_chain": "D"},
    "3c5z": {"antigen_chain": "D"},
    "3rdt": {"antigen_chain": "D"},
    "6dfx": {"antigen_chain": "E"},
    "3c60": {"antigen_chain": "D"},
    "4grl": {"antigen_chain": "D"},
    "6mnn": {"antigen_chain": "D"},
    "6dfs": {"antigen_chain": "D"},
    "4p4k": {"antigen_chain": "B"},
    "4may": {"antigen_chain": "D"},
    "8vd0": {"antigen_chain": "C"},
}

post_fasta_corrections = {
    # 3tf7 has 1 tcrab pair bound together with linker
    "3tf7": {
        "tcr_1_seq": "MGAQSVTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGPQMLLKYYSGDPVVQGVNGFEAEFSKSDSSFHLRKASVHRSDSAVYFCAVSAKGTGSKLSFGKGAKLTVSP",
        "tcr_2_seq": "SEAAVTQSPRNKVTVTGENVTLSCRQTNSHNYMYWYRQDTGHELRLIYYSYGAGNLQIGDVPDGYKATRTTQEDFFLTLESASPSQTSLYFCASSDAPGQLYFGEGSKLTVLELEHHHHHH",
    },
    "8vd0": {"peptide": "GQVELGGGNAVEVCKG"},
    "7q9b": {
        "mhc_chain1": "FFF",
        "mhc_chain2": "GGG",
        "Achain": "III",
        "Bchain": "JJJ",
        "antigen_chain": "HHH",
    },
}

# this is a record of PDB IDs post AF3 cutoff that we find acceptably formatted for our use
post_af3_cutoff_valid_pdbs = [
    "8gom",
    "8vd0",
    "8trr",
    "8wte",
    "8vcy",
    "8es9",
    "8gon",
    "8i5d",
    "8i5c",
    "8vcx",
    "7q99",
    "8eo8",
    "8dnt",
    "8enh",
    "8ye4",
    "8wul",
    "8f5a",
    "8vd2",
    "7q9b",
    "8en8",
    "8pjg",
    "7q9a",
]

# this PDB ID was not found by STCR, but was included in Phil Bradley's 2023 paper
outlier_pdb = pl.DataFrame(
    {
        "pdb": "6l9l",
        "Bchain": "D",
        "Achain": "C",
        "mhc_chain1": "A",
        "mhc_chain2": None,
        "antigen_chain": "B",
        "mhc_class": "I",
        "mhc_1_species": "mouse",
        "mhc_2_species": None,
        "tcr_1_species": None,
        "tcr_2_species": None,
        "mhc_1_chain": "heavy",
        "mhc_2_chain": None,
        "cognate": True,
        "tcr_1_chain": "alpha",
        "tcr_2_chain": "beta",
        "pdb_date": extract_pdb_date({"pdb": "6l9l"})["pdb_date"],
    }
).with_columns(pl.col("pdb_date").str.to_datetime().alias("pdb_date"))


def refmt_rep_dset(rep):
    rep = rep.rename({"pdbid": "pdb"}).with_columns(
        pl.when(pl.col("mhc_class") == 1)
        .then(pl.lit("I"))
        .otherwise(pl.lit("II"))
        .alias("mhc_class"),
    )

    return rep


def refmt_stcr_dest(stcr):
    stcr = format_pdb_df(stcr)
    stcr = get_pdb_date(stcr)

    return stcr


def refmt_mhc_name(mhc_name_df):
    mhc_name_df_II = (
        mhc_name_df.filter(pl.col("mhc_class") == "II")
        .with_columns(pl.col("mhc").str.split(",").alias("split_parts"))
        .with_columns(
            pl.when(pl.col("split_parts").list.len() == 2)
            .then(
                pl.struct(
                    pl.col("split_parts")
                    .list.get(0, null_on_oob=True)
                    .alias("mhc_1_name"),
                    pl.col("split_parts")
                    .list.get(1, null_on_oob=True)
                    .alias("mhc_2_name"),
                )
            )
            .otherwise(
                pl.struct(
                    pl.lit(None).alias("mhc_1_name"),
                    pl.col("split_parts").list.get(0).alias("mhc_2_name"),
                )
            )
            .alias("mhc_struct")
        )
        .unnest("mhc_struct")
    )

    mhc_name_df_I = (
        mhc_name_df.filter(pl.col("mhc_class") == "I")
        .with_columns(pl.col("mhc").str.split(",").alias("split_parts"))
        .with_columns(
            pl.when(pl.col("split_parts").list.len() == 2)
            .then(
                pl.struct(
                    pl.col("split_parts")
                    .list.get(0, null_on_oob=True)
                    .alias("mhc_1_name"),
                    pl.col("split_parts")
                    .list.get(1, null_on_oob=True)
                    .alias("mhc_2_name"),
                )
            )
            .otherwise(
                pl.struct(
                    pl.lit("B2M").alias("mhc_2_name"),
                    pl.col("split_parts").list.get(0).alias("mhc_1_name"),
                )
            )
            .alias("mhc_struct")
        )
        .unnest("mhc_struct")
    )

    return pl.concat(
        [
            mhc_name_df_II,
            mhc_name_df_I,
        ],
        how="vertical",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--raw_csv_path",
        type=str,
    )
    parser.add_argument(
        "-d",
        "--raw_stcr_path",
        type=str,
    )
    parser.add_argument(
        "-i",
        "--imgt_hla_path",
        type=str,
    )
    parser.add_argument(
        "-op",
        "--output_pmhc_path",
        type=str,
    )
    parser.add_argument(
        "-ot",
        "--output_triad_path",
        type=str,
    )
    args = parser.parse_args()

    # https://www.nature.com/articles/s41586-024-07487-w#data-availability
    cutoff = pl.lit(datetime(2023, 1, 12, tzinfo=timezone.utc))

    schema_overrides = {
        "Gchain": pl.String,
        "Dchain": pl.String,
    }
    null_values = ["NA", "unknown", "NOT"]

    stcr = pl.read_csv(
        args.raw_stcr_path,
        schema_overrides=schema_overrides,
        null_values=null_values,
        separator="\t",
    )

    stcr = refmt_stcr_dest(stcr)

    rep = pl.read_csv(
        args.raw_csv_path,
    )

    rep = refmt_rep_dset(rep)

    # we keep a PDB ID for two reasons:
    # 1. it is in Phil Bradley's 2023 paper
    # 2. it is in PDB post-AF3 cutoff and is useable by our pipeline

    rep_pdb_ls = rep.select("pdb").to_series().to_list()

    keep_pdb = rep_pdb_ls + post_af3_cutoff_valid_pdbs

    keep_df = stcr.filter(pl.col("pdb").is_in(keep_pdb))

    keep_df = pl.concat([keep_df, outlier_pdb])

    keep_df = keep_df.with_columns(
        pl.when(pl.col("pdb").is_in(post_af3_cutoff_valid_pdbs))
        .then(pl.lit(False))
        .otherwise(pl.lit(True))
        .alias("replication")
    )

    keep_df = keep_df.join(
        rep.select(
            [
                "pdb",
                "organism",
            ]
        ),
        on="pdb",
        how="left",
    )

    keep_df = keep_df.with_columns(
        pl.when(pl.col("organism").is_null())
        .then(pl.lit("human"))
        .otherwise(pl.col("organism"))
        .alias("tcr_1_species"),
        pl.when(pl.col("organism").is_null())
        .then(pl.lit("human"))
        .otherwise(pl.col("organism"))
        .alias("tcr_2_species"),
    )

    # use the organism information from Phil Bradley's paper
    # the new PDBs happen to all be human

    keep_df = keep_df.with_columns(
        pl.when(pl.col("organism").is_null())
        .then(pl.lit("human"))
        .otherwise(pl.col("organism"))
        .alias("organism"),
    )

    for pdb_id, correction in pre_fasta_corrections.items():

        for k, v in correction.items():
            keep_df = update_df_from_k_v(
                keep_df,
                "pdb",
                pdb_id,
                k,
                v,
            )

    # now extract seqs from FASTA
    keep_df = format_seqs(keep_df)

    # in some cases, peptides were manually extracted in Phil Bradley's paper
    # we keep these for consistency with the paper
    keep_df = (
        keep_df.join(
            rep.select(
                "pdb",
                "peptide",
            ).rename({"peptide": "peptide_manually_extracted"}),
            how="left",
            on="pdb",
        )
        .with_columns(
            pl.when(pl.col("peptide_manually_extracted").is_not_null())
            .then(pl.col("peptide_manually_extracted"))
            .otherwise(pl.col("peptide"))
            .alias("peptide")
        )
        .drop(
            "peptide_manually_extracted",
        )
    )

    for pdb_id, correction in post_fasta_corrections.items():
        for k, v in correction.items():
            keep_df = update_df_from_k_v(
                keep_df,
                "pdb",
                pdb_id,
                k,
                v,
            )

    rep_mhc_names = refmt_mhc_name(rep.select("pdb", "mhc_class", "mhc"))

    new_mhc_names = serial_apply(
        keep_df.filter(~pl.col("replication")),
        infer_correct_mhc,
        HLASequenceDBConverter(args.imgt_hla_path),
        None,
    ).rename(
        {
            "mhc_1_name": "mhc_1_name_inferred",
            "mhc_2_name": "mhc_2_name_inferred",
        }
    )

    keep_df = keep_df.join(
        rep_mhc_names,
        on="pdb",
        how="left",
    )

    keep_df = keep_df.join(
        new_mhc_names.select(
            [
                "pdb",
                "mhc_1_name_inferred",
                "mhc_2_name_inferred",
            ]
        ),
        on="pdb",
        how="left",
    )

    keep_df = keep_df.with_columns(
        pl.when(pl.col("mhc_1_name").is_null())
        .then(pl.col("mhc_1_name_inferred"))
        .otherwise(pl.col("mhc_1_name"))
        .alias("mhc_1_name"),
        pl.when(pl.col("mhc_2_name").is_null())
        .then(pl.col("mhc_2_name_inferred"))
        .otherwise(pl.col("mhc_2_name"))
        .alias("mhc_2_name"),
    )

    keep_df = keep_df.join(
        rep.select(
            [
                "pdb",
                "cdr_rmsd",
                "cdr_rmsd_af2_full",
                "cdr_rmsd_af2_trim",
            ]
        ),
        on="pdb",
        how="left",
    )

    keep_df = generate_job_name(
        keep_df,
        [
            "peptide",
            "mhc_1_seq",
            "mhc_2_seq",
            "tcr_1_seq",
            "tcr_2_seq",
        ],
    )

    keep_df = serial_apply(
        keep_df,
        extract_tcrdist_cols,
    )

    keep_df = serial_apply(keep_df, remove_peptide_from_chains)

    keep_df.select(
        FORMAT_COLS
        + TCRDIST_COLS
        + [
            "pdb",
            "pdb_date",
            "antigen_chain",
            "mhc_chain1",
            "mhc_chain2",
            "Achain",
            "Bchain",
            "organism",
            "cdr_rmsd",
            "cdr_rmsd_af2_full",
            "cdr_rmsd_af2_trim",
        ]
    ).write_parquet(
        args.output_triad_path,
    )

    pdb_antigen = keep_df.select(FORMAT_ANTIGEN_COLS).unique()
    pdb_antigen = generate_job_name(
        pdb_antigen,
        ["peptide", "mhc_1_seq", "mhc_2_seq"],
    )

    pdb_antigen.select(["job_name"] + FORMAT_ANTIGEN_COLS).write_parquet(
        args.output_pmhc_path,
    )
