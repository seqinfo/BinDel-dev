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
    states = 2
    model: hmm.GaussianHMM = hmm.GaussianHMM(n_components=states, covariance_type="full")
    model.startprob_ = np.array([1 / 4000, 1 - 1 / 4000])

    model.transmat_ = np.array([[0.7, 0.3],
                                [0.3, 0.7]])

    model.means_ = np.array([[-0.25, 1.3], [0, 0.5]])
    model.covars_ = np.tile(np.identity(2), (states, 1, 1))

    df[["HMM"]] = model.predict(df[["ratio", "Mann_Whitney"]].to_numpy())
    return df


pd.read_table(os.path.abspath(args.in_file)) \
    .groupby(["sample", "focus"]) \
    .apply(apply_hmm) \
    .to_csv(os.path.abspath(args.out), sep="\t", index=False)
