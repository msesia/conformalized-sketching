import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
import copy
from functools import lru_cache
from tqdm import tqdm
import sys
import mmh3

from collections import defaultdict, OrderedDict
from scipy.stats.mstats import mquantiles
from scipy.special import comb
from scipy.special import loggamma
from scipy import optimize
from scipy import stats
from sklearn import mixture

from cms.data import WordStream, StreamFile, DP
from cms.utils import sort_dict

def dict_to_list(d):
    v = list()
    for x in d.keys():
        v.extend([x]*d[x])
    return v

def lower_bound_from_cdf(pdf, confidence, randomize=False):
    cdf = np.cumsum(pdf[::-1])
    ll = len(cdf) -1 - np.min(np.where(cdf>=confidence)[0])
    idx_0 = len(cdf) -1 - ll
    cdf_0 = cdf[idx_0]
    if randomize:
        idx_1 = idx_0 - 1
        if idx_1<0:
            cdf_1 = 0
        else:
            cdf_1 = cdf[idx_1]
        p_excess = cdf_0 - confidence
        p_randomize = p_excess / (cdf_0-cdf_1)
        randomize = np.random.rand() <= p_randomize
    if randomize:
        ll = ll+1
    return ll

def plot_1dgmm(means, covs, weights, color="black", xlim=None, label=None, plot_mixture=False):
    K = len(means)
    means = np.array(means).reshape(K,1)
    covs = np.array(covs).reshape(K,1)
    weights = np.array(weights).reshape(K,1)

    x_min = np.min(means-3*covs)
    x_max = np.max(means+3*covs)

    # create necessary things to plot
    x_axis = np.arange(x_min, x_max, 0.1)
    y_axes = [stats.norm.pdf(x_axis, float(means[k]), np.sqrt(float(covs[k])))*weights[k] for k in range(K)]
    y_axes = np.array(y_axes)

    #plt.hist(x, density=True, color='black', bins=np.arange(-100, 100, 1))
    for k in range(K):
        plt.plot(x_axis, y_axes[k], lw=3, ls='dashed', c=color, alpha=0.5, label=label)

    if plot_mixture:
        plt.plot(x_axis, np.mean(y_axes,0), c=color, label=label)

    plt.xlim(x_min, x_max)
    if xlim is not None:
        plt.axvline(x=xlim[0])
        plt.axvline(x=xlim[1])

    plt.xlabel(r"X", fontsize=20)
    plt.ylabel(r"Density", fontsize=20)
    plt.legend()

def compute_mode(x, bins=20):
    counts, bins = np.histogram(x, bins=bins)
    max_bin = np.argmax(counts)
    return bins[max_bin:max_bin+2].mean()


def _choice(probs):
    x = np.random.rand()
    cum = 0
    for i,p in enumerate(probs):
        cum += p
        if x < cum:
            break
    return i

class CMS:
    def __init__(self, d, w, seed=2021, conservative=False):
        self.d = d # Number of hash functions
        self.w = w # Width
        self.seed = seed
        self.hash_functions = [self.__generate_hash_function(self.seed+i) for i in range(self.d)]
        self.count = np.zeros((self.d, self.w), dtype='int32')
        self.true_count = defaultdict(lambda: 0)
        self.conservative = conservative

    def reset(self):
        self.count = np.zeros((self.d, self.w), dtype='int32')
        self.true_count = defaultdict(lambda: 0)

    def __generate_hash_function(self, seed=2021):
        """
        Returns a hash function from a family of pairwise-independent hash
        functions

        """
        return lambda x: mmh3.hash(str(x), seed, signed=False) % self.w

    def heavy_hitters_true(self, gamma):
        n = np.sum(self.count[0])
        cutoff = np.floor(gamma * n)
        heavy_hitters = defaultdict(lambda: 0)
        for x in self.true_count.keys():
            if self.true_count[x] >= cutoff:
                heavy_hitters[x] = self.true_count[x]
        return heavy_hitters

    def heavy_hitters_classical(self, gamma):
        n = np.sum(self.count[0])
        cutoff = np.floor(gamma * n)
        heavy_hitters = defaultdict(lambda: 0)
        for x in self.true_count.keys():
            cx = self.estimate_count(x)
            if cx >= cutoff:
                heavy_hitters[x] = cx
        return heavy_hitters


    def classical_error(self, delta):
        n = np.sum(self.count[0])
        epsilon = np.exp(1)/self.w
        error = np.ceil(n * epsilon).astype(int)
        return error

    def apply_hash(self, x):
        columns = np.zeros((self.d,))
        for row in range(self.d):
            columns[row] = self.hash_functions[row](x)
        return columns.astype(int)

    def update_count(self, x, n=1):
        self.true_count[x] += 1
        columns = self.apply_hash(x)
        c_hat = self.estimate_count(x)

        for row in range(self.d):
            if self.conservative:
                current_count = self.count[row, columns[row]]
                new_count = np.maximum(current_count, c_hat + n)
                self.count[row, columns[row]] = new_count
            else:
                self.count[row, columns[row]] += n

        return columns

    def estimate_count(self, x):
        value = float("inf")
        columns = self.apply_hash(x)
        for row in range(self.d):
            value = np.minimum(value, self.count[row, columns[row]])
        return value.astype(int)

    def lower_bound(self, x, confidence):
        error = self.classical_error(1.0-confidence)
        upper = self.estimate_count(x)
        lower = np.maximum(0, upper - error).astype(int)
        return lower, None

    def print(self):
        od = OrderedDict(sorted(self.true_count.items(), key=lambda x:x[1], reverse=True))

        results = []
        for x in od:
            y_true = od[x]
            y_hat = self.estimate_count(x)
            results.append({'x':x, 'count':y_true, 'upper': y_hat})

        return pd.DataFrame(results)

class BayesianDP:
    def __init__(self, cms, alpha=None, sigma=None, tau=None):
        self.cms = copy.deepcopy(cms)
        self.C = self.cms.count
        self.alpha = alpha

        if hasattr(alpha, "__len__"):
            alpha_list = alpha
        else:
            alpha_list = np.array([alpha])
        self.param_list = alpha_list

    def posterior(self, x, param=None):
        if param is None:
            param = self.alpha

        alpha = param

        N = self.C.shape[0]
        J = self.C.shape[1]
        m = np.sum(self.C[0])

        def _log_pmf_c_k(J, c, k_max):
            k = np.arange(k_max+1)
            out = np.zeros((k_max+1,))
            if np.isinf(alpha):
                out = -np.inf * np.ones((k_max+1,))
                out = 0
            else:
                out = np.log(alpha/J)
                out += loggamma(c + 1.0) - loggamma(c - k + 1.0)
                out += loggamma(c - k + alpha/J) - loggamma(c + alpha/J + 1.0)
            return out

        def _posterior(c_v):
            N = self.C.shape[0]
            K = np.min(c_v)

            log_v = np.zeros((K+1,))
            for n in range(N):
                log_v += _log_pmf_c_k(J, c_v[n], K)

            # Compute the part that was missing from Cai's paper
            if N>1:
                 tmp = _log_pmf_c_k(1, m, K)
                 log_v += (1.0-N) * tmp

            prob = np.exp(log_v - np.max(log_v))
            prob /= np.sum(prob)

            return prob


        columns = self.cms.apply_hash(x)
        c_v = [self.C[row,columns[row]] for row in range(self.C.shape[0])]
        prob = _posterior(c_v)

        return prob

    def lower_bound(self, x, confidence, param=None, randomize=False):
        if param is None:
            param = self.alpha
        pdf = self.posterior(x, param=param)
        ll = lower_bound_from_cdf(pdf, confidence, randomize=randomize)
        return ll, pdf

    def _neg_log_likelihood(self, alpha):
        """Compute negative log-likelihood for given α"""
        N = self.C.shape[0]
        J = self.C.shape[1]
        M = np.sum(self.C,1)[0]

        def prop_log_V_n(alpha, n):
            log_v = loggamma(self.C[n,:] + alpha/J) - loggamma(self.C[n,:]+1) - loggamma(alpha/J)
            log_v = np.sum(log_v)
            return log_v

        def prop_log_V(alpha):
            log_v = 0
            for n in range(N):
                log_v += loggamma(M+1) + loggamma(alpha) - loggamma(M+alpha) + prop_log_V_n(alpha,n)
            return log_v

        ll = prop_log_V(alpha)
        return -ll

    def empirical_bayes(self):
        """Estimate α via Empirical Bayes"""
        J = self.C.shape[1]
        opt_sol = optimize.minimize_scalar(self._neg_log_likelihood, method="Bounded", bounds=(0.0001,100*J))
        self.alpha = opt_sol.x
        return self.alpha

class BayesianCMS:
    def __init__(self, stream, cms, model="DP", alpha=None, sigma=None, tau=None, posterior="mcmc"):
        if alpha is not None:
            assert (alpha>0)
        if sigma is not None:
            assert (sigma>=0) and (sigma<1)
        if tau is not None:
            assert (tau>=0) and (tau<=1)

        self.stream = stream
        self.cms = copy.deepcopy(cms)
        self.alpha = alpha
        self.sigma = sigma
        self.tau = tau
        self.model_name = model
        self.posterior = posterior

    def run(self, n, n_test, confidence=0.9, seed=2021):

        np.random.seed(seed)

        self.cms.reset()
        ## Process stream
        true_frequency = defaultdict(lambda: 0)
        self.stream.reset()
        print("Processing training data....")
        sys.stdout.flush()
        for i in tqdm(range(n)):
            x = self.stream.sample()
            self.cms.update_count(x)

        # Initialize Bayesian model
        if self.model_name=="DP":
            model = BayesianDP(self.cms, alpha=self.alpha)

            if self.alpha is None:
                alpha_hat = model.empirical_bayes()
                print("Empirical Bayes estimated parameter: {:.3f}".format(alpha_hat))

        else:
            print("Error! Unknown model.")
            pdb.set_trace()

        # Evaluate
        print("Evaluating on test data....")
        sys.stdout.flush()
        np.random.seed(seed)
        results = []
        for i in tqdm(range(n_test)):
            x = self.stream.sample()
            y = self.cms.true_count[x]
            upper = self.cms.estimate_count(x)

            # Compute lower bound using exact posterior
            lower, _ = model.lower_bound(x, confidence, randomize=True)
            posterior = model.posterior(x)

            # Estimate the true count
            post_mean = np.sum(np.arange(len(posterior))*posterior)
            post_median = np.min(np.where(np.cumsum(posterior)>0.5)[0])
            post_mode = np.argmax(posterior)

            # Return results
            results.append({'method':'Bayesian',
                            'x':x, 'count':y, 'upper': upper, 'lower':lower,
                            'mean':post_median, 'median':post_median, 'mode':post_mode,
                            'seen':False})

        results = pd.DataFrame(results)
        results["tracking"] = False
        results = results.sort_values(by=['count'], ascending=False)
        return results

class ClassicalCMS:
    def __init__(self, stream, cms):
        self.stream = stream
        self.cms = copy.deepcopy(cms)

    def run(self, n, n_test, confidence=0.9, seed=None):

        self.cms.reset()
        ## Process stream
        true_frequency = defaultdict(lambda: 0)
        self.stream.reset()
        for i in tqdm(range(n)):
            x = self.stream.sample()
            self.cms.update_count(x)

        # Evaluate
        print("Evaluating on test data....")
        sys.stdout.flush()
        results = []
        for i in tqdm(range(n_test)):
            x = self.stream.sample()
            y = self.cms.true_count[x]
            upper = self.cms.estimate_count(x)
            lower, _ = self.cms.lower_bound(x, confidence)
            mean = (upper+lower)/2
            results.append({'method': 'Classical', 'x':x, 'count':y, 'upper': upper, 'lower':lower,
                            'mean':mean, 'median':mean, 'mode':mean,
                            'seen':False})

        results = pd.DataFrame(results)
        results = results.sort_values(by=['count'], ascending=False)
        return results
