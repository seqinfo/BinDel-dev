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

# TODO: GMMHMM ? https://hmmlearn.readthedocs.io/en/latest/api.html#hmmlearn.hmm.GMMHMM
# TODO: andmeteaduse meetodid
# TODO: Gaussi segumudel hoopis?
def apply_hmm(df: pd.DataFrame):
    states: int = 3
    model: hmm.GaussianHMM = hmm.GaussianHMM(n_components=states, covariance_type="full")
    model.startprob_ = np.array([0, 1, 0])
    model.transmat_ = np.array(states * [states * [1 / states]])
    model.means_ = np.array([[-0.7, 1.2], [0, 0.45], [0.7, 1.2]])
    model.covars_ = np.tile(np.identity(2), (states, 1, 1))

    df[["HMM"]] = model.predict(df[["ratio", "Mann_Whitney"]].to_numpy())
    return df


pd.read_table(os.path.abspath(args.in_file)) \
    .groupby(["sample", "focus"]) \
    .apply(apply_hmm) \
    .to_csv(os.path.abspath(args.out), sep="\t", index=False)
