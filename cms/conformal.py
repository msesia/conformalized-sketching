import sys
import numpy as np
import pandas as pd
import pdb
import time

from collections import defaultdict, OrderedDict, Counter
from scipy.stats.mstats import mquantiles
from scipy.stats import mannwhitneyu
from scipy.stats import norm
from scipy.stats import t as student_t
import copy

from cms.cms import BayesianCMS, BayesianDP
from cms.cqr import QR, QRScores
from cms.utils import sum_dict, dictToList, listToDict
from cms.chr import HistogramAccumulator

from sklearn.ensemble import IsolationForest
from sklearn.svm import OneClassSVM
from sklearn.neighbors import LocalOutlierFactor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold

from tqdm import tqdm

import matplotlib.pyplot as plt

from functools import lru_cache

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
    if (a_set & b_set):
        return True
    else:
        return False

def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

class ClassicalScores:
    def __init__(self, cms, method="constant"):
        self.cms = copy.deepcopy(cms)
        self.method = method
        delta_t = 1
        self.n = np.sum(self.cms.count[0])

    def compute_score(self, x, y):
        "This score measures by how much we need to decrease the upper bound to obtain a valid lower bound"
        upper = self.cms.estimate_count(x)
        score = upper-y
        return score, 0

    def predict_interval(self, x, tau, tau_u):
        upper = self.cms.estimate_count(x)
        lower = np.maximum(0, upper - tau).astype(int)
        return lower, upper

class ClassicalScoresTwoSided:
    def __init__(self, cms, method="constant"):
        self.cms = copy.deepcopy(cms)
        self.method = method
        delta_t = 1
        self.n = np.sum(self.cms.count[0])

    def compute_score(self, x, y):
        "This score measures by how much we need to decrease the upper bound to obtain a valid lower bound"
        upper = self.cms.estimate_count(x)
        lower = upper
        score_low = lower-y
        score_upp = y-upper
        return score_low, score_upp

    def predict_interval(self, x, tau_l, tau_u):
        upper = self.cms.estimate_count(x)
        lower = upper
        lower = np.maximum(0, lower - tau_l).astype(int)
        upper = np.maximum(lower, np.minimum(upper, upper + tau_u)).astype(int)
        return lower, upper

class BayesianScores:
    def __init__(self, model, confidence):
        self.model = model
        self.confidence = confidence
        self.t_seq = np.linspace(0, 1, 100)

    def _lower_bound(self, cdfi, t):
        lower = np.where(cdfi>=t)[0]
        if len(lower) > 0:
            lower = np.max(lower)
        else:
            lower = 0
        return lower

    def _compute_sequence(self, x):
        posterior = self.model.posterior(x)
        cdfi = np.cumsum(posterior[::-1])
        t_seq = self.t_seq.reshape((len(self.t_seq),1))
        A = cdfi >= t_seq
        lower_seq = np.sum(A,1).astype(int)
        lower_seq[lower_seq<0] = 0
        return lower_seq

    @lru_cache()
    def compute_score(self, x, y):
        lower = self._compute_sequence(x)
        idx_below = np.where(lower <= y)[0]
        if len(idx_below)>0:
            score = np.min(idx_below)
        else:
            score = len(lower)-1

        return score, 0

    def predict_interval(self, x, t, t_u):
        upper = self.cms.estimate_count(x)
        posterior = self.model.posterior(x)
        cdfi = np.cumsum(posterior[::-1])
        tau = self.t_seq[t]
        idx_above = np.where(cdfi>=tau-1e-6)[0]
        lower = len(posterior) - np.min(idx_above) - 1
        return lower, upper

class BayesianScoresTwoSided:
    def __init__(self, model, confidence):
        self.model = model
        self.confidence = confidence
        self.t_seq = np.linspace(0, 1, 100)

    def _lower_bound(self, cdfi, t):
        lower = np.where(cdfi>=t)[0]
        if len(lower) > 0:
            lower = np.max(lower)
        else:
            lower = 0
        return lower

    @lru_cache()
    def compute_score(self, x, y):
        pdf = self.model.posterior(x)
        pdf = pdf.reshape((1,len(pdf)))
        breaks = np.arange(pdf.shape[1])
        CHR = HistogramAccumulator(pdf, breaks, self.confidence, delta_alpha=0.01)
        score = CHR.calibrate_intervals(y)
        return score, 0

    def predict_interval(self, x, t, t_u):
        pdf = self.model.posterior(x)
        pdf = pdf.reshape((1,len(pdf)))
        breaks = np.arange(pdf.shape[1])
        CHR = HistogramAccumulator(pdf, breaks, self.confidence, delta_alpha=0.01)
        S = CHR.predict_intervals(1.0-t)
        return S.flatten().astype(int)

class BootstrapScores:
    def __init__(self, cms, alpha):
        self.cms = copy.deepcopy(cms)
        self.alpha = alpha

    @lru_cache()
    def estimate_noise_dist(self, x, n_mc = 1000):
        r,w = self.cms.count.shape
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
    def compute_score(self, x, y):
        "This score measures by how much we need to decrease the upper bound to obtain a valid lower bound"
        upper_max = self.cms.estimate_count(x)
        noise = self.estimate_noise_dist(x)
        delta_max = mquantiles(noise, 1.0-self.alpha)[0]
        lower = np.maximum(0, upper_max - delta_max)
        score_low = lower-y
        return score_low, 0

    @lru_cache()
    def predict_interval(self, x, tau_l, tau_u):
        upper_max = self.cms.estimate_count(x)
        noise = self.estimate_noise_dist(x)
        delta_max = mquantiles(noise, 1.0-self.alpha)[0]
        lower = np.maximum(0, upper_max - delta_max)
        upper = upper_max
        lower = np.maximum(0, lower - tau_l).astype(int)
        return lower, upper


class BootstrapScoresTwoSided:
    def __init__(self, cms, alpha, n_mc = 1000):
        self.cms = copy.deepcopy(cms)
        self.n_mc = n_mc
        self.alpha = alpha

    @lru_cache()
    def estimate_noise_dist(self, x, n_mc = 1000):
        r,w = self.cms.count.shape
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
    def compute_score(self, x, y):
        "This score measures by how much we need to decrease the upper bound to obtain a valid lower bound"
        upper_max = self.cms.estimate_count(x)
        noise = self.estimate_noise_dist(x)
        delta_min = mquantiles(noise, self.alpha/2)[0]
        delta_max = mquantiles(noise, 1.0-self.alpha/2)[0]
        lower = np.maximum(0, upper_max - delta_max)
        upper = np.maximum(0, upper_max - delta_min)
        score_low = lower-y
        score_upp = y-upper
        return score_low, score_upp

    @lru_cache()
    def predict_interval(self, x, tau_l, tau_u):
        upper_max = self.cms.estimate_count(x)
        noise = self.estimate_noise_dist(x)
        delta_min = mquantiles(noise, self.alpha/2)[0]
        delta_max = mquantiles(noise, 1.0-self.alpha/2)[0]
        lower = np.maximum(0, upper_max - delta_max)
        upper = np.maximum(0, upper_max - delta_min)
        lower = np.maximum(0, lower - tau_l).astype(int)
        upper = np.maximum(lower, np.minimum(upper, upper + tau_u)).astype(int)
        return lower, upper


class BootstrapScoresTwoSidedCHR:
    def __init__(self, cms, confidence, n_mc = 1000):
        self.cms = copy.deepcopy(cms)
        self.alpha = confidence
        self.n_mc = n_mc
        self.t_seq = np.linspace(0, 1, 100)

    @lru_cache()
    def _estimate_noise_dist(self, x):
        # CMS parameters
        r,w = self.cms.count.shape

        noise = np.zeros((self.n_mc,))
        i = 0
        while i < self.n_mc:
            I = self.cms.apply_hash(x)
            J = np.random.choice(w, r)
            if not common_member(I,J):
                c_J = np.array([self.cms.count[k, J[k]] for k in range(r)])
                noise[i] = np.min(c_J)
                i = i +1
        return noise

    @lru_cache()
    def estimate_median(self, x):
        upper = self.cms.estimate_count(x)
        noise = self._estimate_noise_dist(x).astype(int)
        vals = np.maximum(upper - noise, 0)
        return np.median(vals)

    @lru_cache()
    def estimate_quantiles(self, x):
        upper = self.cms.estimate_count(x)
        noise = self._estimate_noise_dist(x).astype(int)
        vals = np.maximum(upper - noise, 0)
        return mquantiles(vals, [self.alpha, 1.0-self.alpha])

    @lru_cache()
    def compute_score(self, x, y):
        lower, upper = self.estimate_quantiles(x)
        score = np.maximum(lower-y, y-upper)
        return score, 0

    @lru_cache()
    def predict_interval(self, x, t_l, t_u):
        upper_max = self.cms.estimate_count(x)
        lower, upper = self.estimate_quantiles(x)
        lower = lower - t_l
        upper = upper + t_u
        lower = np.maximum(0, lower)
        upper = np.minimum(upper_max, upper)
        S = np.array([lower,upper])
        return S.astype(int)


class ConformalCMS:
    def __init__(self, stream, cms, n_track, prop_train=0.5, n_bins=1, scorer_type="Bayesian-DP", two_sided=False, unique=1):
        self.stream = stream
        self.cms = cms
        self.max_track = n_track
        self.prop_train = prop_train
        self.n_bins = n_bins
        self.scorer_type = scorer_type
        self.unique = unique
        self.two_sided = two_sided

    @lru_cache()
    def _predict_interval(self, x, scorer=None, t_hat_low=None, t_hat_upp=None):
        lower_warmup = self.cms_warmup.true_count[x]
        if scorer is not None:
            lower, upper = scorer.predict_interval(x, t_hat_low, t_hat_upp)
            if hasattr(lower, "__len__"):
                lower = lower[0]
                upper = upper[0]
            lower = np.maximum(0, lower)
        else:
            lower = 0
            upper = self.cms.estimate_count(x)

        return lower + lower_warmup, upper + lower_warmup

    def run(self, n, n_test, confidence=0.9, seed=2021, heavy_hitters_gamma=0.01, shift=0):
        n_bins = self.n_bins
        scorer_type = self.scorer_type

        print("Running conformal method with n = {:d}...".format(n))
        sys.stdout.flush()

        ## Warmup
        print("Warm-up iterations (max track: {:d})...".format(self.max_track))
        sys.stdout.flush()
        freq_track = defaultdict(lambda: 0)
        data_track = defaultdict(lambda: 0)
        self.cms_warmup = copy.deepcopy(self.cms)
        #n_track = 0
        n_track = self.max_track
        i_range=tqdm(range(self.max_track))
        for i in i_range:
            x = self.stream.sample()
            self.cms_warmup.update_count(x)
            freq_track[x] = 0
            data_track[x] += 1

        n1 = n - n_track
        print("Main iterations: {:d}...".format(n1))
        sys.stdout.flush()
        # Process stream
        for i in tqdm(range(n1)):
            x = self.stream.sample()
            self.cms.update_count(x)

            # Check whether this object is being tracked
            if x in freq_track.keys():
                freq_track[x] += 1

        if n1 > 0:
            print("Calibrating lower bound....")
            sys.stdout.flush()

            J = self.cms.w

            if scorer_type == "Bayesian-DP":
                model = BayesianDP(self.cms)
                alpha_hat = model.empirical_bayes()
                if self.two_sided:
                    scorer = BayesianScoresTwoSided(model, 1.0-confidence)
                else:
                    scorer = BayesianScores(model, 1.0-confidence)

            elif scorer_type == "Bootstrap":
                if self.two_sided:
                    scorer = BootstrapScoresTwoSided(self.cms, 1.0-confidence)
                else:
                    scorer = BootstrapScores(self.cms, 1.0-confidence)

            elif scorer_type == "Constant":
                if self.two_sided:
                    scorer = ClassicalScoresTwoSided(self.cms, method="constant")
                else:
                    scorer = ClassicalScores(self.cms, method="constant")

            elif scorer_type == "Proportional":
                scorer = ClassicalScores(self.cms, method="proportional")

            elif scorer_type == "Adaptive":
                scorer = QRScores(self.cms, confidence, seed=seed, two_sided=self.two_sided)
                #scorer = HRScores(self.cms, confidence, seed=seed)

                # Split the tracking data into training/calibration sets
                if self.unique > 1:
                    Y_track = np.array(list(freq_track.values()))
                    object_track = np.array(list(freq_track.keys()))
                    upper_track_list = np.array([self.cms.estimate_count(x) for x in freq_track.keys()])
                    freq_track_keys = np.array(list(freq_track.keys()))
                else:
                    Y_track = np.concatenate([[freq_track[x]]*data_track[x] for x in freq_track.keys()])
                    object_track = np.concatenate([[x]*data_track[x] for x in freq_track.keys()])
                    upper_track_list = np.concatenate([[self.cms.estimate_count(x)]*data_track[x] for x in freq_track.keys()])
                    freq_track_keys = np.concatenate([[x]*data_track[x] for x in freq_track.keys()])

                upper_track = upper_track_list.reshape((len(upper_track_list),1))
                upper_train, upper_calib, Y_train, Y_calib, _, freq_track_keys_calib = train_test_split(upper_track, Y_track, freq_track_keys,
                                                                                                        test_size=1-self.prop_train, random_state=2022)
                # Train the predictive model
                scorer.train(upper_train, Y_train, X_calib=upper_calib, Y_calib=Y_calib)
                # Repackage the calibration data
                freq_track = defaultdict(lambda: 0)
                data_track = defaultdict(lambda: 0)
                for i in range(len(freq_track_keys_calib)):
                    x = freq_track_keys_calib[i]
                    freq_track[x] = Y_calib[i]
                    data_track[x] = np.sum(freq_track_keys_calib==x)

            else:
                raise("Error! Unknown type of scores.")

            if self.unique > 1:
                scores_cal_tmp = np.concatenate([[scorer.compute_score(x, freq_track[x])]*data_track[x] for x in tqdm(freq_track.keys())])
                y_cal_tmp = np.concatenate([[freq_track[x]]*data_track[x] for x in freq_track.keys()])
                x_cal_tmp = np.concatenate([[x]*data_track[x] for x in freq_track.keys()])
                # Split calibration data points into subsets
                G = int(np.maximum(2,np.floor(len(y_cal_tmp)/self.unique)))
                scores_cal = np.zeros((G,2))
                y_cal = np.zeros((G,))
                kf = KFold(n_splits=G, random_state=2022, shuffle=True)
                idx_fold_list = [idx for _, idx in kf.split(x_cal_tmp)]
                for g in range(G):
                    idx_fold = idx_fold_list[g]
                    _, idx_unique = np.unique(x_cal_tmp[idx_fold], return_index=True)
                    idx_pick = np.random.choice(idx_unique)
                    scores_cal[g] = scores_cal_tmp[idx_fold][idx_pick]
                    y_cal[g] = y_cal_tmp[idx_fold][idx_pick]

            else:
                # Calculate scores
                scores_cal = np.concatenate([[scorer.compute_score(x, freq_track[x])]*data_track[x] for x in tqdm(freq_track.keys())])
                y_cal = np.concatenate([[freq_track[x]]*data_track[x] for x in freq_track.keys()])

            # Calibrate the conformity scores (for bin-conditional coverage)
            n_bins_max = int(np.maximum(1, np.floor(len(y_cal)/100)))
            n_bins = np.minimum(n_bins, n_bins_max)
            y_bins, y_bin_cutoffs = pd.qcut(y_cal, n_bins, duplicates="drop", labels=False, retbins=True, precision=0)
            print("Cutoffs for {:d} bins:".format(n_bins))
            print(y_bin_cutoffs)

            calibrated_scores_bins_low = [None]*n_bins
            calibrated_scores_bins_upp = [None]*n_bins
            for k in range(n_bins):
                idx_bin = np.where(y_bins==k)[0]
                n_bin = len(idx_bin)
                alpha = 1.0 - confidence
                if len(idx_bin) > 0:
                    if self.two_sided:
                        level_adjusted = (1.0-alpha/2)*(1.0+1.0/float(n_bin))
                    else:
                        level_adjusted = (1.0-alpha)*(1.0+1.0/float(n_bin))

                    calibrated_scores_bins_low[k] = mquantiles(scores_cal[idx_bin,0], prob=level_adjusted)[0]
                    calibrated_scores_bins_upp[k] = mquantiles(scores_cal[idx_bin,1], prob=level_adjusted)[0]
                    if not self.two_sided:
                        calibrated_scores_bins_low[k] = np.ceil(calibrated_scores_bins_low[k]).astype(int)
                        calibrated_scores_bins_upp[k] = np.ceil(calibrated_scores_bins_upp[k]).astype(int)
                else:
                    calibrated_scores_bins_low[k] = 0
                    calibrated_scores_bins_upp[k] = 0
            calibrated_score_low = np.max(calibrated_scores_bins_low)
            calibrated_score_upp = np.max(calibrated_scores_bins_upp)
            print("Calibrated scores (low):")
            print(calibrated_scores_bins_low)
            print("Calibrated scores (upp):")
            print(calibrated_scores_bins_upp)
            print("Calibrated score (final):")
            print([calibrated_score_low, calibrated_score_upp])

            # Combine warm-up and regular cms
            self.cms.count = self.cms.count + self.cms_warmup.count
            self.cms.true_count = sum_dict(self.cms.true_count, self.cms_warmup.true_count)
        else:
            scorer = None
            calibrated_score_low = None
            calibrated_score_upp = None
            # Combine warm-up and regular cms
            self.cms.count = self.cms_warmup.count
            self.cms.true_count = self.cms_warmup.true_count


        # Evaluate
        print("Evaluating on test data....")
        sys.stdout.flush()
        np.random.seed(seed)
        results = []
        for i in tqdm(range(n_test)):
            x = self.stream.sample()
            if shift > 0:
                if np.random.rand() < shift:
                    x = x + np.random.rand()
            y = self.cms.true_count[x]
            lower, upper = self._predict_interval(x, scorer=scorer, t_hat_low=calibrated_score_low, t_hat_upp=calibrated_score_upp)
            tracking = x in freq_track

            # Estimation
            if confidence==0.5:
                est_median = lower
            else:
                est_median = (upper+lower)/2

            results.append({'method':'Conformal-'+scorer_type,
                            'x':x, 'count':y, 'upper': upper, 'lower':lower,
                            'mean':est_median, 'median':est_median, 'mode':est_median,
                            'seen':tracking})

        results = pd.DataFrame(results)
        results = results.sort_values(by=['count'], ascending=False)
        return results
