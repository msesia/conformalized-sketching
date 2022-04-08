import sys
import numpy as np
import pandas as pd
import pdb

from collections import defaultdict, OrderedDict, Counter
from scipy.stats.mstats import mquantiles
import copy

from cms.cms import BayesianCMS, BayesianDP
from cms.cqr import QR, QRScores, HRScores
from cms.utils import sum_dict, dictToList, listToDict

from sklearn.model_selection import train_test_split

from tqdm import tqdm

import matplotlib.pyplot as plt


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
        return score

    def predict_lower(self, x, tau):
        upper = self.cms.estimate_count(x)
        lower = np.maximum(0, upper - tau).astype(int)
        return lower

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

    def compute_score(self, x, y):
        lower = self._compute_sequence(x)
        idx_below = np.where(lower <= y)[0]
        if len(idx_below)>0:
            score = np.min(idx_below)
        else:
            score = len(lower)-1

        return score

    def predict_lower(self, x, t):
        posterior = self.model.posterior(x)
        cdfi = np.cumsum(posterior[::-1])
        tau = self.t_seq[t]
        idx_above = np.where(cdfi>=tau-1e-6)[0]
        lower = len(posterior) - np.min(idx_above) - 1
        return lower

class ConformalCMS:
    def __init__(self, stream, cms, n_track, prop_train=0.5, n_bins=1, scorer_type="Bayesian-DP", unique=False):
        self.stream = stream
        self.cms = cms
        self.max_track = n_track
        self.prop_train = prop_train
        self.n_bins = n_bins
        self.scorer_type = scorer_type
        self.unique = unique

    # def _predict(self, x, alpha, scorer, scores_cal, upper_cal, n_bins=1, t_hat=None):
    #     if t_hat is None:
    #         upper = self.cms.estimate_count(x)
    #         upper_aug = np.concatenate([upper_cal, [upper]])
    #         upper_bins = pd.qcut(upper_aug, n_bins)
    #         upper_bins_cal = upper_bins[:-1]
    #         upper_bin = upper_bins[-1]
    #         idx_bin = np.where(upper_bins_cal==upper_bin)[0]
    #         n_cal = len(scores_cal)
    #         if n_cal > 0:
    #             level_adjusted = (1.0-alpha)*(1.0+1.0/float(n_cal))
    #             calibrated_score = (mquantiles(scores_cal[idx_bin], prob=level_adjusted)[0]).astype(int)
    #             t_hat = calibrated_score
    #         else:
    #             t_hat = -1
    #     lower = scorer.predict_lower(x, t_hat)
    #     lower = np.maximum(0, lower)
    #     lower_warmup = self.cms_warmup.true_count[x]
    #     return lower + lower_warmup

    def _predict(self, x, scorer=None, t_hat=None):
        lower_warmup = self.cms_warmup.true_count[x]
        if scorer is not None:
            lower = scorer.predict_lower(x, t_hat)
            if hasattr(lower, "__len__"):
                lower = lower[0]
            lower = np.maximum(0, lower)
        else:
            lower = 0
        return lower + lower_warmup

    def run(self, n, n_test, confidence=0.9, seed=2021, heavy_hitters_gamma=0.01):
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
                scorer = BayesianScores(model, 1.0-confidence)

            elif scorer_type == "Constant":
                scorer = ClassicalScores(self.cms, method="constant")

            elif scorer_type == "Proportional":
                scorer = ClassicalScores(self.cms, method="proportional")

            elif scorer_type == "Adaptive":
                scorer = QRScores(self.cms, confidence, seed=seed)
                #scorer = HRScores(self.cms, confidence, seed=seed)

                # Split the tracking data into training/calibration sets
                if self.unique:
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

            if self.unique:
                # This approach is not being used.
                # Compute dictionary of conformity scores
                num_scores = sum(freq_track.values())
                scores_unique = defaultdict(lambda: 0)
                scores_dict = defaultdict(lambda: 0)
                for x in tqdm(freq_track.keys()):
                    y = freq_track[x]
                    score = scorer.compute_score(x, y)
                    scores_dict[score] += data_track[x]
                    scores_unique[x] = score
                scores_cal = np.array(list(scores_unique.values()))
                y_cal = np.array(list(freq_track.values()))
            else:
                scores_cal = np.concatenate([[scorer.compute_score(x, freq_track[x])]*data_track[x] for x in freq_track.keys()])
                y_cal = np.concatenate([[freq_track[x]]*data_track[x] for x in freq_track.keys()])

            # Calibrate the conformity scores (for bin-conditional coverage)
            n_bins_max = int(np.maximum(1, np.floor(len(y_cal)/100)))
            n_bins = np.minimum(n_bins, n_bins_max)
            y_bins, y_bin_cutoffs = pd.qcut(y_cal, n_bins, duplicates="drop", labels=False, retbins=True, precision=0)
            print("Cutoffs for {:d} bins:".format(n_bins))
            print(y_bin_cutoffs)

            calibrated_scores_bins = [None]*n_bins
            for k in range(n_bins):
                idx_bin = np.where(y_bins==k)[0]
                n_bin = len(idx_bin)
                alpha = 1.0 - confidence
                if len(idx_bin) > 0:
                    level_adjusted = (1.0-alpha)*(1.0+1.0/float(n_bin))
                    calibrated_scores_bins[k] = np.ceil(mquantiles(scores_cal[idx_bin], prob=level_adjusted)[0]).astype(int)
                else:
                    calibrated_scores_bins[k] = 0
            calibrated_score = np.max(calibrated_scores_bins)
            print("Calibrated scores:")
            print(calibrated_scores_bins)
            print("Calibrated score (final):")
            print(calibrated_score)

            # Combine warm-up and regular cms
            self.cms.count = self.cms.count + self.cms_warmup.count
            self.cms.true_count = sum_dict(self.cms.true_count, self.cms_warmup.true_count)
        else:
            scorer = None
            calibrated_score = None
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
            y = self.cms.true_count[x]
            upper = self.cms.estimate_count(x)
            lower = self._predict(x, scorer=scorer, t_hat=calibrated_score)
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
