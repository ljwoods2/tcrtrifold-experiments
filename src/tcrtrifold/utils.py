import polars as pl
import hashlib


def hash_sequence(seq: str, hash_type: str = "md5") -> str:
    """
    Hash a TCR sequence using the specified hash function.

    Args:
        tcr_seq (str): The TCR sequence string.
        hash_type (str): The hash function to use ('md5', 'sha1', 'sha256', etc.)

    Returns:
        str: The hexadecimal digest of the hashed sequence.
    """
    # Select the hash function
    if hash_type.lower() == "md5":
        h = hashlib.md5()
    elif hash_type.lower() == "sha1":
        h = hashlib.sha1()
    elif hash_type.lower() == "sha256":
        h = hashlib.sha256()
    else:
        raise ValueError("Unsupported hash type")

    # Encode the sequence and compute the hash
    h.update(seq.encode("utf-8"))
    return h.hexdigest()


def generate_job_name(df, cols, name="job_name"):
    df = df.with_columns(
        pl.concat_str(
            pl.concat_str(
                [
                    *[pl.col(colname) for colname in cols],
                ],
                ignore_nulls=True,
            )
            .map_elements(
                lambda x: hash_sequence(x, "md5"), return_dtype=pl.String
            )
            .alias(name),
        )
    )
    return df


def update_df_from_k_v(
    df,
    primary_key_colname,
    primary_key,
    k,
    v,
):
    df = pl.concat(
        [
            df.filter(pl.col(primary_key_colname) == primary_key).with_columns(
                pl.lit(v).alias(k)
            ),
            df.filter(pl.col(primary_key_colname) != primary_key),
        ],
        how="vertical_relaxed",
    )
    return df
