#!/usr/bin/env python
"""


"""
from tcrtrifold.utils import generate_job_name
from tcrtrifold.pdb import *
from tcr_format_parsers.common.TriadUtils import (
    generate_job_name,
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
}

post_fasta_corrections = {
    # 3tf7 has 1 tcrab pair bound together with linker
    "3tf7": {
        "tcr_1_seq": "MGAQSVTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGPQMLLKYYSGDPVVQGVNGFEAEFSKSDSSFHLRKASVHRSDSAVYFCAVSAKGTGSKLSFGKGAKLTVSP",
        "tcr_2_seq": "SEAAVTQSPRNKVTVTGENVTLSCRQTNSHNYMYWYRQDTGHELRLIYYSYGAGNLQIGDVPDGYKATRTTQEDFFLTLESASPSQTSLYFCASSDAPGQLYFGEGSKLTVLELEHHHHHH",
    }
}


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

    all_pdb_stcr = format_pdb_df(
        pl.read_csv(
            args.raw_stcr_path,
            schema_overrides=schema_overrides,
            null_values=null_values,
            separator="\t",
        )
    )

    all_pdb_stcr = serial_apply(
        all_pdb_stcr,
        get_pdb_date,
    )

    rep = (
        pl.read_csv(
            "raw/table_S1_structure_benchmark_complexes.csv",
        )
        .rename({"pdbid": "pdb"})
        .with_columns(
            pl.when(pl.col("mhc_class") == 1)
            .then(pl.lit("I"))
            .otherwise(pl.lit("II"))
            .alias("mhc_class"),
        )
    )
