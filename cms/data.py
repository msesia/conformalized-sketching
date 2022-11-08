import numpy as np
from numpy.random import RandomState
import matplotlib.pyplot as plt
from scipy import stats
from tqdm import tqdm
import csv
import sys

from sklearn.datasets import fetch_20newsgroups
from sklearn.feature_extraction.text import CountVectorizer

from collections import defaultdict, OrderedDict

import nltk
from nltk.util import ngrams
from nltk.tokenize import sent_tokenize
tokenizer = nltk.RegexpTokenizer(r"\w+")
from english_words import english_words_set

from Bio import SeqIO
import screed

import pdb

def random_choice(probs, rng=None):
    if rng is None:
        x = np.random.rand()
    else:
        x = rng.random()
    cum = 0
    for i,p in enumerate(probs):
        cum += p
        if x < cum:
            break
    return i

class Zipf:
    def __init__(self, alpha, seed=2021):
        assert (alpha>0)
        self.alpha = alpha
        self.seed = seed
        self.reset()

    def reset(self, seed=None):
        if seed is not None:
            self.seed = seed

        # Initialize random number generator
        self.rng = np.random.Generator(np.random.PCG64(self.seed))

    def sample(self, n=1, verbose=False, method=None, store=None):
        # Sample sequentially
        data = self.rng.zipf(self.alpha, n)
        if n==1:
            data=data[0]
        return data

class PYP:
    """ Pitman-Yor Process"""
    def __init__(self, alpha, sigma, seed=2021):
        assert (sigma>=0) and (sigma<1)
        self.alpha = alpha
        self.sigma = sigma
        self.seed = seed
        self.reset()

    def reset(self, seed=None):
        if seed is not None:
            self.seed = seed
        # Initialize data container
        self.data = defaultdict(lambda: 0)

        # Initialize random number generator
        self.rng = np.random.Generator(np.random.PCG64(self.seed))
        self.P0 = self.rng.standard_normal

    def _prob_vec(self, counts):
        # Compute probability of observing a new label or a previously seen label
        m = np.sum(counts)
        k = len(counts)
        # Compute probability of observing a new label
        if m>0:
            p_unseen = (self.alpha+k*self.sigma) / (self.alpha + m)
        else:
            p_unseen = 1.0
        # Compute probability of observing a previously seen label
        if m>0:
            p_seen = (counts-self.sigma) / (self.alpha + m)
        else:
            p_seen = 0 * (counts-self.sigma)
        out = np.maximum(0, np.concatenate([[p_unseen], p_seen]))
        if np.abs(1-np.sum(out))>1e-3:
            pdb.set_trace()
        return out/np.sum(out)

    def _sample_step(self, data):
        counts = np.fromiter(data.values(), dtype=int)
        prob = self._prob_vec(counts)
        #j = np.random.choice(len(prob), 1, p=prob)[0]
        j = random_choice(prob, rng=self.rng)
        if j == 0:
            x = np.round(self.P0(), 6)
        else:
            x = list(data.keys())[j-1]
        return x

    def sample_k(self, n):
        k = 0
        for m in np.arange(n):
            # Compute probability of observing a new label
            if m>0:
                p_unseen = self.sigma * k / m
            else:
                p_unseen = 1
            if np.random.uniform() < p_unseen:
                k += 1
        return k

    def sample(self, n=1, verbose=False, store=True, method=None):
        # Initialize empty data container
        data = defaultdict(lambda: 0)

        # Sample sequentially
        for i in range(n):
            x = self._sample_step(self.data)
            data[x] = data[x]+1
            if store:
                self.data[x] = self.data[x] + 1

        if n==1:
            data = x

        return data

class SP:
    """ Stable process """
    def __init__(self, sigma, seed=2021):
        assert (sigma>=0) and (sigma<1)
        self.sigma = sigma
        self.seed = seed
        self.reset()

    def reset(self, seed=None):
        if seed is not None:
            self.seed = seed
        # Initialize data container
        self.data = defaultdict(lambda: 0)

        # Initialize random number generator
        self.rng = np.random.Generator(np.random.PCG64(self.seed))
        self.P0 = self.rng.standard_normal

    def _prob_vec(self, counts):
        # Compute probability of observing a new label or a previously seen label
        m = np.sum(counts)
        k = len(counts)
        # Compute probability of observing a new label
        if m>0:
            p_unseen = self.sigma * k / m
        else:
            p_unseen = 1.0
        # Compute probability of observing a previously seen label
        if m>0:
            p_seen = (counts-self.sigma) / m
        else:
            p_seen = 0 * (counts-self.sigma)
        out = np.maximum(0, np.concatenate([[p_unseen], p_seen]))
        if np.abs(1-np.sum(out))>1e-3:
            pdb.set_trace()
        return out/np.sum(out)

    def _sample_step(self, data):
        counts = np.fromiter(data.values(), dtype=int)
        prob = self._prob_vec(counts)
        #j = np.random.choice(len(prob), 1, p=prob)[0]
        j = random_choice(prob, rng=self.rng)
        if j == 0:
            x = np.round(self.P0(), 6)
        else:
            x = list(data.keys())[j-1]
        return x

    def sample_k(self, n):
        k = 0
        for m in np.arange(n):
            # Compute probability of observing a new label
            if m>0:
                p_unseen = self.sigma * k / m
            else:
                p_unseen = 1
            if np.random.uniform() < p_unseen:
                k += 1
        return k

    def sample(self, n=1, verbose=False, store=True, method=None):
        # Initialize empty data container
        data = defaultdict(lambda: 0)

        # Sample sequentially
        for i in range(n):
            x = self._sample_step(self.data)
            data[x] = data[x]+1
            if store:
                self.data[x] = self.data[x] + 1

        if n==1:
            data = x

        return data

class DP:
    def __init__(self, alpha, seed=2021):
        assert (alpha>0)
        self.alpha = alpha
        self.seed = seed
        self.reset()

    def reset(self, seed=None):
        if seed is not None:
            self.seed = seed

        # Initialize data container
        self.data = defaultdict(lambda: 0)

        # Initialize random number generator
        self.rng = np.random.Generator(np.random.PCG64(self.seed))
        self.P0 = self.rng.standard_normal

    def _prob_vec(self, counts):
        # Compute probability of observing a new label or a previously seen label
        m = np.sum(counts)
        k = len(counts)
        # Add place holder for new label
        counts = np.concatenate([[self.alpha], counts])
        # Compute probabilities for new and old labels
        prob = counts / (self.alpha + m)
        assert np.abs(np.sum(prob) - 1) < 1e-3
        return prob/np.sum(prob)

    def _sample_step(self, data):
        counts = np.fromiter(data.values(), dtype=int)
        prob = self._prob_vec(counts)
        #j = np.random.choice(len(prob), 1, p=prob)[0]
        j = random_choice(prob, rng=self.rng)
        if j == 0:
            x = np.round(self.P0(), 6)
        else:
            x = list(data.keys())[j-1]
        return x

    def sample(self, n=1, verbose=False, store=True, method=None):
        # Initialize empty data container
        data = defaultdict(lambda: 0)

        # Sample sequentially
        for i in range(n):
            x = self._sample_step(self.data)
            data[x] = data[x]+1
            if store:
                self.data[x] = self.data[x] + 1
        if n==1:
            data = x
        return data

def extract_ngrams(s, n, words):
    s = s.lower()
    sentences = sent_tokenize(s)
    tokens_raw = [tokenizer.tokenize(t) for t in sentences]
    tokens = [[w for w in t if w in words] for t in tokens_raw]
    n_grams = [ngrams(t, n) for t in tokens]
    data_ngrams = np.concatenate([[ ' '.join(grams) for grams in n_grams[i]] for i in range(len(n_grams))])
    return data_ngrams

class WordStream:
    def __init__(self, n_docs=100, n_grams=2, seed=2021, filename_out=None):
        file_list = nltk.corpus.gutenberg.fileids()
        n_docs = np.minimum(n_docs, len(file_list))
        self.data = []
        for k in tqdm(range(n_docs)):
            filename = file_list[k]
            data_raw = nltk.corpus.gutenberg.raw(filename)
            data_ngrams = extract_ngrams(data_raw, n_grams, english_words_set)
            self.data.append(data_ngrams)
        self.data = np.concatenate(self.data)

        if filename_out is not None:
            with open(filename_out,'w') as data_file:
                wr = csv.writer(data_file, dialect='excel')
                wr.writerow(self.data)
                print("List of {:d} tuples written to: {:s}".format(len(self.data), filename_out))
                sys.stdout.flush()

        self.counter = 0
        self.prng = RandomState(seed)

    def set_seed(self, seed):
        self.prng = RandomState(seed)

    def sample(self, n=1, verbose=None):
        if n==1:
            idx = self.prng.randint(0, len(self.data))
            words = self.data[idx]
        else:
            words = [None]*n
            for i in range(n):
                idx = self.prng.randint(0, len(self.data))
                words[i] = self.data[idx]
        return words

    def reset(self):
        return None

class DNAStream:
    def __init__(self, k_size=10, seed=2021, seq_max=100, filename_out=None):
        def build_kmers(sequence, ksize):
            kmers = []
            n_kmers = len(sequence) - ksize + 1

            for i in range(n_kmers):
                kmer = sequence[i:i + ksize]
                kmers.append(kmer)

            return kmers

        def read_kmers_from_file(filename_in, filename_out, ksize, seq_max=100):
            all_kmers = []
            with open(filename_out,'w') as out_file:
                wr = csv.writer(out_file, dialect='excel')
                with screed.open(filename_in) as seqfile:
                    i = 0
                    for record in tqdm(seqfile):
                        if i >= seq_max:
                            break
                        else:
                            i+=1
                        kmers = build_kmers(record.sequence, ksize)
                        wr.writerow(kmers)

        #self.data = read_kmers_from_file("data/clear_covid_ca.fasta", k_size)
        if filename_out is not None:
            read_kmers_from_file("data/clear_covid_ca.fasta", filename_out, k_size, seq_max=seq_max)
            print("Finished writing file {:s}.".format(filename_out))
            sys.stdout.flush()

        self.counter = 0
        self.prng = RandomState(seed)

    def set_seed(self, seed):
        self.prng = RandomState(seed)

    def sample(self, n=1, verbose=None):
        if n==1:
            idx = self.prng.randint(0, len(self.data))
            words = self.data[idx]
        else:
            words = [None]*n
            for i in range(n):
                idx = self.prng.randint(0, len(self.data))
                words[i] = self.data[idx]
        return words

    def reset(self):
        return None


class StreamFile(WordStream):
    def __init__(self, filename, seed=2021):
        with open(filename,'r') as data_file:
            rd = csv.reader(data_file, dialect='excel')
            self.data = list(rd)
            self.data = np.concatenate(self.data)

        print("Loaded {:d} data points.".format(len(self.data)))
        self.counter = 0
        self.prng = RandomState(seed)
