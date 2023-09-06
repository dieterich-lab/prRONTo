#!/usr/bin/env python

import click
import pandas as pd
import numpy as np
from sklearn.neighbors import LocalOutlierFactor


@click.command()
@click.option("-f", "--feature", multiple=True, required=True, help="Columns to use as score")
@click.option("-o", "--output", required=True, help="Output file name")
@click.option("-n", "--neighbors", default=20, help="Number of neighbors for LOF calculation")
@click.option("-c", "--contamination", default=0.002, help="Contamination value for LOF calculation")
@click.option("-m", "--filter-median", help="Filter outlier beyond median")
@click.argument("file", type=click.Path(exists=True))
def cli(scores, output, neighbors, contamination, filter_median, file):
  df = pd.read_csv(file, sep = "\t")
  lof_params = {
          "neighbors": neighbors,
          "contamination": contamination,}
  new_df = add_lof_score(df, list(scores), filter_median, lof_params)
  new_df.to_csv(output, sep="\t", index=False)


def add_lof_score(df, cols, filter_median = False, lof_params = None):
  lof_params = lof_params if  lof_params else {}
  new_col_lof = "_".join(cols) + "_lof_score"
  new_col_outlier = "_".join(cols) + "_lof_outlier"

  lof = LocalOutlierFactor(contamination=lof_params.get("contamination", "auto"),
                           n_neighbors=lof_params.get("neighbors", 20))
  ohat = lof.fit_predict(df[cols])
  scores = lof.negative_outlier_factor_
  df[new_col_lof] = scores
  df[new_col_outlier] = ohat
  if filter_median:
      smaller = df[cols].median() > df[cols]
      # reduce rows
      reduced_rows= [all(row) for row in smaller.itertuples(index=False)]
      df[reduced_rows, new_col_outlier] = 0

  return df


if __name__ == "__main__":
    cli()
