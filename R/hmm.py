import os
from argparse import ArgumentParser

import numpy as np
import pandas as pd
from hmmlearn import hmm

np.random.seed(0)

parser: ArgumentParser = ArgumentParser()

parser.add_argument("-i", "--input", dest="in_file", help="input file")
parser.add_argument("-o", "--output", dest="out", help="output file name")

args = parser.parse_args()


def apply_hmm(df: pd.DataFrame):
    model: hmm.GaussianHMM = hmm.GaussianHMM(n_components=3, covariance_type="full")
    model.startprob_ = np.array([0.25, 0.5, 0.25])

    model.transmat_ = np.array([[1 / 3, 1 / 3, 1 / 3],
                                [1 / 3, 1 / 3, 1 / 3],
                                [1 / 3, 1 / 3, 1 / 3]])

    model.means_ = np.array([[-0.4, 0.7], [0, 0.25], [0.4, 0.7]])
    model.covars_ = np.tile(np.identity(2), (3, 1, 1))

    df[["HMM"]] = model.predict(df[["ratio", "Mann_Whitney"]].to_numpy())
    return df


pd.read_table(os.path.abspath(args.in_file)) \
    .groupby(["sample", "focus"]) \
    .apply(apply_hmm) \
    .to_csv(os.path.abspath(args.out), sep="\t", index=False)
