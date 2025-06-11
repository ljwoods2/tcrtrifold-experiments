import requests
from Bio import SeqIO
from io import StringIO
import polars as pl
from tcr_format_parsers.common.MHCCodeConverter import (
    HLASequenceDBConverter,
    H2SequenceDictConverter,
)
from tcr_format_parsers.common.TCRUtils import standardize_tcr
import warnings
from pathlib import Path
import MDAnalysis as mda

SEQ_STRUCT = pl.Struct(
    {
        "peptide_seq": pl.String,
        "mhc_1_seq": pl.String,
        "mhc_2_seq": pl.String,
        "tcr_1_seq": pl.String,
        "tcr_2_seq": pl.String,
    }
)


def format_pdb_df(df):
    df = df.with_columns(
        pl.when(pl.col("mhc_type") == "MH1")
        .then(pl.lit("I"))
        .when(pl.col("mhc_type") == "MH2")
        .then(pl.lit("II"))
        .otherwise(None)
        .alias("mhc_class"),
    )

    # df = df.filter(
    #     (pl.col("mhc_chain1").is_not_null())
    #     & (pl.col("mhc_chain2").is_not_null())
    # )

    df = df.group_by("pdb").agg(
        pl.col("Bchain").drop_nulls().first(),
        pl.col("Achain").drop_nulls().first(),
        pl.col("mhc_chain1").drop_nulls().first(),
        pl.col("mhc_chain2").drop_nulls().first(),
        pl.col("antigen_chain").drop_nulls().first(),
        pl.col("mhc_class").drop_nulls().first(),
        pl.col("mhc_chain1_organism").drop_nulls().first().alias("mhc_1_species"),
        pl.col("mhc_chain2_organism").drop_nulls().first().alias("mhc_2_species"),
        pl.col("alpha_organism").drop_nulls().first().alias("tcr_1_species"),
        pl.col("beta_organism").drop_nulls().first().alias("tcr_2_species"),
    )

    df = df.with_columns(
        pl.when(pl.col("mhc_1_species") == "homo sapiens")
        .then(pl.lit("human"))
        .when(pl.col("mhc_1_species") == "mus musculus")
        .then(pl.lit("mouse"))
        .otherwise(None)
        .alias("mhc_1_species"),
        pl.when(pl.col("mhc_2_species") == "homo sapiens")
        .then(pl.lit("human"))
        .when(pl.col("mhc_2_species") == "mus musculus")
        .then(pl.lit("mouse"))
        .otherwise(None)
        .alias("mhc_2_species"),
        pl.when(pl.col("tcr_1_species") == "homo sapiens")
        .then(pl.lit("human"))
        .when(pl.col("tcr_1_species") == "mus musculus")
        .then(pl.lit("mouse"))
        .otherwise(None)
        .alias("tcr_1_species"),
        pl.when(pl.col("tcr_2_species") == "homo sapiens")
        .then(pl.lit("human"))
        .when(pl.col("tcr_2_species") == "mus musculus")
        .then(pl.lit("mouse"))
        .otherwise(None)
        .alias("tcr_2_species"),
    )

    df = df.with_columns(
        pl.when(pl.col("mhc_class") == "II")
        .then(pl.lit("alpha"))
        .otherwise(pl.lit("heavy"))
        .alias("mhc_1_chain"),
        pl.when(pl.col("mhc_class") == "II")
        .then(pl.lit("beta"))
        .otherwise(pl.lit("light"))
        .alias("mhc_2_chain"),
        pl.lit(True).alias("cognate"),
        pl.lit("alpha").alias("tcr_1_chain"),
        pl.lit("beta").alias("tcr_2_chain"),
    )

    df = df.with_columns(
        pl.col("antigen_chain")
        .str.split("|")
        .list.first()
        .str.strip_chars()
        .alias("antigen_chain")
    )

    return df


def get_pdb_date(row):
    r = requests.get("https://data.rcsb.org/rest/v1/core/entry/" + row["pdb"])
    r.raise_for_status()
    new_row = row.copy()
    new_row["pdb_date"] = r.json()["rcsb_accession_info"][
        "initial_release_date"
    ]

    return pl.DataFrame(new_row).with_columns(
        pl.col("pdb_date").str.to_datetime().alias("pdb_date")
    )


def parse_chain(chain):
    if "[" in chain:
        
        return chain.split("[auth ")[1][0]
        # if can have multi-letter chains
        # return chain.split("[auth ")[1].split("]")[0]
    else:
        return chain.replace(" ", "")


def parse_fasta_description(description):
    chain_token = description.split("|")[1]

    if chain_token.startswith("Chain "):
        return list(parse_chain(chain_token.split("Chain ")[1]))
    else:
        chains = chain_token.split("Chains ")[1].split(",")
        chain_list = [parse_chain(chain) for chain in chains]

        return chain_list


def get_fasta_seq(
    pdb_id,
    antigen_chain_id,
    mhc_chain1_id,
    mhc_chain2_id,
    Achain_id,
    Bchain_id,
):
    r = requests.get("https://www.rcsb.org/fasta/entry/" + pdb_id)

    r.raise_for_status()

    fasta_sequences = SeqIO.parse(StringIO(r.text), "fasta")

    seq_dict = {}
    for fasta in fasta_sequences:
        chains = parse_fasta_description(fasta.description)
        for chain in chains:
            seq_dict[chain] = str(fasta.seq)

    return {
        "peptide_seq": seq_dict[antigen_chain_id] if antigen_chain_id is not None else "",
        "mhc_1_seq": seq_dict[mhc_chain1_id] if mhc_chain1_id is not None else "",
        "mhc_2_seq": seq_dict[mhc_chain2_id] if mhc_chain2_id is not None else "",
        "tcr_1_seq": seq_dict[Achain_id] if Achain_id is not None else "",
        "tcr_2_seq": seq_dict[Bchain_id] if Bchain_id is not None else "",
    }



def format_seqs(df, skip_peptide=False):
    df = df.with_columns(
        pl.struct(
            pl.col("pdb"),
            pl.col("Bchain"),
            pl.col("Achain"),
            pl.col("antigen_chain"),
            pl.col("mhc_chain1"),
            pl.col("mhc_chain2"),
        )
        .map_elements(
            lambda x: get_fasta_seq(
                x["pdb"],
                x["antigen_chain"],
                x["mhc_chain1"],
                x["mhc_chain2"],
                x["Achain"],
                x["Bchain"],
            ),
            return_dtype=SEQ_STRUCT,
            skip_nulls=False,
        )
        .alias("chain_seqs"),
    ).unnest("chain_seqs").with_columns(
        pl.when(pl.col("peptide_seq") == "")
        .then(pl.lit(None))
        .otherwise(pl.col("peptide_seq"))
        .alias("peptide_seq"),
        pl.when(pl.col("mhc_1_seq") == "")
        .then(pl.lit(None))
        .otherwise(pl.col("mhc_1_seq"))
        .alias("mhc_1_seq"),
        pl.when(pl.col("mhc_2_seq") == "")
        .then(pl.lit(None))
        .otherwise(pl.col("mhc_2_seq"))
        .alias("mhc_2_seq"),
        pl.when(pl.col("tcr_1_seq") == "")
        .then(pl.lit(None))
        .otherwise(pl.col("tcr_1_seq"))
        .alias("tcr_1_seq"),
        pl.when(pl.col("tcr_2_seq") == "")
        .then(pl.lit(None))
        .otherwise(pl.col("tcr_2_seq"))
        .alias("tcr_2_seq"),
    )

    return df


def remove_peptide_from_chains(row):
    new_row = row.copy()

    if row["mhc_1_seq"] is not None and row["peptide"] in row["mhc_1_seq"]:
        warnings.warn(f"Peptide found in MHC 1 sequence for PDB {row['pdb']} at position {row['mhc_1_seq'].index(row['peptide'])}")
        index_of_peptide = row["mhc_1_seq"].index(row["peptide"])
        new_row["mhc_1_seq"] = new_row["mhc_1_seq"][
            index_of_peptide + len(row["peptide"]) :
        ]
    if row["mhc_2_seq"] is not None and row["peptide"] in row["mhc_2_seq"]:
        warnings.warn(f"Peptide found in MHC 2 sequence for PDB {row['pdb']} at position {row['mhc_2_seq'].index(row['peptide'])}")
        index_of_peptide = row["mhc_2_seq"].index(row["peptide"])
        new_row["mhc_2_seq"] = new_row["mhc_2_seq"][
            index_of_peptide + len(row["peptide"]) :
        ]
    if row["tcr_1_seq"] is not None and row["peptide"] in row["tcr_1_seq"]:
        warnings.warn(f"Peptide found in TCR 1 sequence for PDB {row['pdb']} at position {row['tcr_1_seq'].index(row['peptide'])}")
        index_of_peptide = row["tcr_1_seq"].index(row["peptide"])
        new_row["tcr_1_seq"] = new_row["tcr_1_seq"][
            index_of_peptide + len(row["peptide"]) :
        ]
    if row["tcr_2_seq"] is not None and row["peptide"] in row["tcr_2_seq"]:
        warnings.warn(f"Peptide found in TCR 2 sequence for PDB {row['pdb']} at position {row['tcr_2_seq'].index(row['peptide'])}")
        index_of_peptide = row["tcr_2_seq"].index(row["peptide"])
        new_row["tcr_2_seq"] = new_row["tcr_2_seq"][
            index_of_peptide + len(row["peptide"]) :
        ]
    return pl.DataFrame(new_row)


def infer_correct_mhc(row, human_conv, mouse_conv):
    mhc1 = row["mhc_1_seq"]
    mhc2 = row["mhc_2_seq"]

    if row["organism"] == "human":
        mhc_1_inf = human_conv.get_mhc_allele(
            mhc1, chain=row["mhc_1_chain"], top_only=True
        )
    else:
        mhc_1_inf = mouse_conv.get_mhc_allele(
            mhc1, chain=row["mhc_1_chain"], top_only=True
        )

    if row["organism"] == "human":
        if row["mhc_class"] == "I" and row['mhc_2_seq'] is None:
            mhc_2_inf = {
                "mhc_2_match_seq": None,
                "mhc_2_name" : None,
                "mhc_2_match_size": None,
                "mhc_2_match_proportion": None,
                "mhc_2_status": None,
                "mhc_2_maxres": None,
            }
        else:
            mhc_2_inf = human_conv.get_mhc_allele(
                mhc2, chain=row["mhc_2_chain"], top_only=True
            )
    else:
        mhc_2_inf = mouse_conv.get_mhc_allele(
            mhc2, chain=row["mhc_2_chain"], top_only=True
        )

    new_row = row.copy()

    new_row["mhc_1_match_seq"] = mhc_1_inf["seq"]
    new_row["mhc_1_name"] = mhc_1_inf["name"]
    new_row["mhc_1_match_size"] = mhc_1_inf["match_size"]
    new_row["mhc_1_match_proportion"] = (
        (mhc_1_inf["match_size"] / len(mhc1))
        if mhc_1_inf["match_size"] is not None
        else None
    )
    new_row["mhc_1_status"] = mhc_1_inf["sequence_status"]
    new_row["mhc_1_name_maxres"] = mhc_1_inf["max_resolution_name"]

    new_row["mhc_2_match_seq"] = mhc_2_inf["seq"]
    new_row["mhc_2_name"] = mhc_2_inf["name"]
    new_row["mhc_2_match_size"] = mhc_2_inf["match_size"]
    new_row["mhc_2_match_proportion"] = (
        (mhc_2_inf["match_size"] / len(mhc2))
        if mhc_2_inf["match_size"] is not None
        else None
    )
    new_row["mhc_2_status"] = mhc_2_inf["sequence_status"]
    new_row["mhc_2_maxres"] = mhc_2_inf["max_resolution_name"]
    return pl.DataFrame(new_row)


def download_pdb(row, path):

    r = requests.get(f"https://files.rcsb.org/download/{row["pdb"]}.pdb")
    suffix = ".pdb"
    try:
        r.raise_for_status()
    except Exception as e:
        r = requests.get(f"https://files.rcsb.org/download/{row["pdb"]}.cif")
        suffix = ".cif"
        r.raise_for_status()
    with open(path / (row["pdb"] + suffix), "wb") as f:
        f.write(r.content)
    return pl.DataFrame(row)

def get_true_mda_universe(pdb_id, root_path):
    # Favor PDB since it doesn't have multiple residue with same ID issue
    if (root_path / (pdb_id + ".pdb")).exists():
        suffix = ".pdb"

    else:
        suffix = ".cif"

    return mda.Universe((root_path / (pdb_id + suffix)).as_posix())

