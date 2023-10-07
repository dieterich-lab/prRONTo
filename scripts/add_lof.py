#!/usr/bin/env python

import click
import pandas as pd
import numpy as np
from sklearn.neighbors import LocalOutlierFactor


@click.command()
@click.option("-f", "--feature", type=str, multiple=True, required=True, help="Columns to use as score")
@click.option("-o", "--output", type=str, required=True, help="Output file name")
@click.option("-n", "--neighbors", type=int, default=20, help="Number of neighbors for LOF calculation")
@click.option("-c", "--contamination", type=float, default=0.002, help="Contamination value for LOF calculation")
@click.option("-m", "--filter-median", is_flag=True, show_default=True, default=False, help="Filter outlier beyond median")
@click.argument("file", type=click.Path(exists=True))
def cli(feature, output, neighbors, contamination, filter_median, file):
  df = pd.read_csv(file, sep = "\t")
  lof_params = {
          "neighbors": neighbors,
          "contamination": contamination,}
  new_df = add_lof_score(df, list(feature), filter_median, lof_params)
  new_df.to_csv(output, sep="\t", index=False)


def add_lof_score(df, features, filter_median, lof_params):
  for feature in features:
    new_col_lof = "lof_score_" + feature
    new_col_outlier =  "lof_outlier_" + feature
    cols = [f"feature_{col}" for col in feature.split("_")]

    lof = LocalOutlierFactor(contamination=lof_params["contamination"],
                             n_neighbors=lof_params["neighbors"])
    ohat = lof.fit_predict(df[cols])
    scores = lof.negative_outlier_factor_
    df[new_col_lof] = scores * -1 # make scores positive
    df[new_col_outlier] = -1 # flip class labels: -1 -> inlier

    smaller = df[cols].median() > df[cols]
    # reduce
    smaller = smaller.eval("&".join(smaller))
    sorted_scores = ~smaller[np.argsort(scores)] # sorted scores and boolean if bigger than median
    n_outliers = round(len(scores) * lof_params["contamination"])
    if filter_median:
      # own calculation of outlier to filter by median
      i_outliers = sorted_scores.index[sorted_scores == True]
    else:
      i_outliers = sorted_scores.index
    # FIXME
    df.loc[i_outliers[0:n_outliers], new_col_outlier] = 1

  return df


if __name__ == "__main__":
    cli()
