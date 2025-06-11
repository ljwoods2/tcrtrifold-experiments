def sample_to(antigen_df, neg_df, neg_per_pos):
    """
    For each antigen in antigen_df, select neg_per_pos negatives from neg_df

    Negative selection follows these rules:

    1. If there aren't enough negatives in neg_df, add all of them, and add
        the antigen to the returned missing_df with a needed_negatives column
    2. If there are enough negatives in neg_df, select neg_per_pos negatives.
        Select using the sampling method:
        1. Select a random source antigen
        2. Select a random TCR from that source antigen
        3. Remove that TCR from all source antigens (since they may be shared across source antigens)
    """

    new_neg_df = []
    missing_df = []

    for row in tqdm(
        antigen_df.iter_rows(named=True),
        total=antigen_df.height,
        desc="Processing rows",
    ):
        antigen = pl.DataFrame([row])
        desired_neg = row["TCRdiv_samples"] * neg_per_pos

        # remove source, which is duplicated
        true_neg_df = (
            neg_df.join(antigen, on=FORMAT_ANTIGEN_COLS)
            .select(FORMAT_COLS + TCRDIST_COLS)
            .unique()
        )
        true_neg = true_neg_df.height

        if desired_neg > true_neg:
            needed_neg = desired_neg - true_neg
            # jsut add them all
            new_neg_df.append(true_neg_df)
            missing_df.append(
                antigen.with_columns(
                    pl.lit(needed_neg).alias("needed_negatives")
                )
            )

        elif true_neg > desired_neg:

            candidate_neg = neg_df.join(antigen, on=FORMAT_ANTIGEN_COLS)
            neg_by_source_antigen = candidate_neg.partition_by(
                SOURCE_ANTIGEN_COLS
            )
            neg_partition_meta = [df.height for df in neg_by_source_antigen]
            remaining_neg = list(range(len(neg_by_source_antigen)))

            selected = 0
            sel_tcrs = []

            while selected < desired_neg:

                ant_int = random.choice(remaining_neg)

                # neg_partition_meta[ant_int] -= 1
                selected += 1

                # sample a TCR and strike it out from all source antigen
                sel_tcr = (
                    neg_by_source_antigen[ant_int]
                    .select(FORMAT_TCR_COLS + TCRDIST_COLS)
                    .unique()
                    .sample(n=1, shuffle=True)
                )
                sel_tcrs.append(sel_tcr)

                rmvlist = []
                for i in remaining_neg:
                    # strikeout TCR
                    neg_by_source_antigen[i] = neg_by_source_antigen[i].join(
                        sel_tcr, on=FORMAT_TCR_COLS + TCRDIST_COLS, how="anti"
                    )

                    # update remaining TCRs to select
                    neg_partition_meta[i] = neg_by_source_antigen[i].height

                    if neg_partition_meta[i] == 0:
                        rmvlist.append(i)

                for i in rmvlist:
                    remaining_neg.remove(i)

            sel_tcrs = pl.concat(sel_tcrs)

            noncognate = antigen.join(sel_tcrs, how="cross")
            noncognate = generate_job_name(noncognate)
            noncognate = noncognate.with_columns(
                pl.lit(False).alias("cognate")
            ).select(FORMAT_COLS + TCRDIST_COLS)

            if noncognate.height != noncognate.unique().height:
                print(row)

            new_neg_df.append(noncognate)

        else:
            # exactly the right number of negs
            new_neg_df.append(true_neg_df)

    return pl.concat(new_neg_df, how="vertical_relaxed"), pl.concat(
        missing_df, how="vertical_relaxed"
    )


def sample_supplemental_negatives(missing_neg_antigens, cross_class_negatives):
    """
    Takes the missing_neg_antigens df from sample_to and the cross_class_negatives
    which is a dataframe of all possible negative TCRs from the opposite MHC class.

    For each antigen in missing_neg_antigens, select the number of negatives
    needed from the cross_class_negatives. Use the same sampling method as
    sample_to, but only select from the cross_class_negatives.
    """
    new_neg_df = []
    for row in tqdm(
        missing_neg_antigens.iter_rows(named=True),
        total=missing_neg_antigens.height,
        desc="Processing rows",
    ):
        antigen = pl.DataFrame([row])

        needed_neg = row["needed_negatives"]

        candidate_neg = cross_class_negatives.join(
            antigen, on=FORMAT_ANTIGEN_COLS
        )
        neg_by_source_antigen = candidate_neg.partition_by(SOURCE_ANTIGEN_COLS)
        neg_partition_meta = [df.height for df in neg_by_source_antigen]
        remaining_neg = list(range(len(neg_by_source_antigen)))

        if needed_neg > candidate_neg.height:
            raise ValueError("Not enough negatives")

        selected = 0
        sel_tcrs = []

        while selected < needed_neg:

            ant_int = random.choice(remaining_neg)

            # neg_partition_meta[ant_int] -= 1
            selected += 1

            # sample a TCR and strike it out from all source antigen
            sel_tcr = (
                neg_by_source_antigen[ant_int]
                .select(FORMAT_TCR_COLS + TCRDIST_COLS)
                .unique()
                .sample(n=1, shuffle=True)
            )
            sel_tcrs.append(sel_tcr)

            rmvlist = []
            for i in remaining_neg:
                # strikeout TCR
                neg_by_source_antigen[i] = neg_by_source_antigen[i].join(
                    sel_tcr, on=FORMAT_TCR_COLS + TCRDIST_COLS, how="anti"
                )

                # update remaining TCRs to select
                neg_partition_meta[i] = neg_by_source_antigen[i].height

                if neg_partition_meta[i] == 0:
                    rmvlist.append(i)

            for i in rmvlist:
                remaining_neg.remove(i)

        sel_tcrs = pl.concat(sel_tcrs)

        noncognate = antigen.join(sel_tcrs, how="cross")
        noncognate = generate_job_name(noncognate)
        noncognate = noncognate.with_columns(
            pl.lit(False).alias("cognate")
        ).select(FORMAT_COLS + TCRDIST_COLS)

        if noncognate.height != noncognate.unique().height:
            print(row)

        new_neg_df.append(noncognate)

    return pl.concat(new_neg_df, how="vertical_relaxed")
