import numpy as np
from statsmodels.stats.multitest import multipletests
import sys
import pandas as pd

from collections import defaultdict, OrderedDict, Counter
from scipy.stats.mstats import mquantiles
from scipy.stats import mannwhitneyu

from cms.cms import BayesianCMS, BayesianDP
from cms.cqr import QR, QRScores
from cms.conformal import ClassicalScores, BayesianScores
from cms.utils import sum_dict, dictToList, listToDict
from cms.chr import HistogramAccumulator

import copy

from sklearn.ensemble import IsolationForest
from sklearn.svm import OneClassSVM
from sklearn.neighbors import LocalOutlierFactor
from sklearn.model_selection import train_test_split

from tqdm import tqdm

import matplotlib.pyplot as plt

import pdb

from functools import lru_cache

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
    if (a_set & b_set):
        return True
    else:
        return False


class BootstrapCMS:
    def __init__(self, stream, cms, two_sided=False):
        self.stream = stream
        self.cms = cms
        self.two_sided = two_sided

    def run(self, n, n_test, confidence=0.9, seed=2021, shift=0):
        print("Running bootstrap method with n = {:d}...".format(n))
        sys.stdout.flush()

        # Process stream
        for i in tqdm(range(n)):
            x = self.stream.sample()
            self.cms.update_count(x)

        # CMS parameters
        r,w = self.cms.count.shape

        @lru_cache()
        def estimate_noise_dist(x, n_mc = 10000):
            noise = np.zeros((n_mc,))
            i = 0
            while i < n_mc:
                I = self.cms.apply_hash(x)
                J = np.random.choice(w, r)
                if not common_member(I,J):
                    c_J = np.array([self.cms.count[k, J[k]] for k in range(r)])
                    noise[i] = np.min(c_J)
                    i = i +1
            return noise

        @lru_cache()
        def compute_pdf(x):
            upper = self.cms.estimate_count(x)
            noise = estimate_noise_dist(x).astype(int)
            vals = np.maximum(upper - noise, 0)
            pdf = np.zeros((upper+1,))
            vals_bins = np.bincount(vals)
            pdf[np.arange(len(vals_bins))] = vals_bins
            pdf = pdf + 0.01 / len(pdf)
            pdf = pdf / np.sum(pdf)
            pdf = pdf.reshape((1,len(pdf)))
            return pdf

        # Evaluate
        print("Evaluating on test data....")
        sys.stdout.flush()
        np.random.seed(seed)
        results = []
        noise = None
        for i in tqdm(range(n_test)):
            x = self.stream.sample()
            y = self.cms.true_count[x]

            if noise is None:
                noise = estimate_noise_dist(x)

            upper_max = self.cms.estimate_count(x)
            if self.two_sided:
                alpha = 1.0 - confidence
                delta_min = mquantiles(noise, alpha/2)[0]
                delta_max = mquantiles(noise, 1.0-alpha/2)[0]
                lower = np.maximum(0, upper_max - delta_max)
                upper = np.maximum(0, upper_max - delta_min)
                    
            else:
                delta = mquantiles(noise, confidence)[0]
                lower = np.maximum(0, upper_max - delta)
                upper = upper_max

            # Estimation
            delta_median = mquantiles(noise, 0.5)[0]
            est_median = np.maximum(0, upper_max - delta_median)

            results.append({'method':'Bootstrap',
                            'x':x, 'count':y, 'upper': upper, 'lower':lower,
                            'mean':est_median, 'median':est_median, 'mode':est_median,
                            'seen':0})

        results = pd.DataFrame(results)
        results = results.sort_values(by=['count'], ascending=False)
        return results
