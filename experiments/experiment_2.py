import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pdb

import os, sys
sys.path.append("..")

from cms.data import WordStream, StreamFile, DP, SP, PYP, Zipf
from cms.cms import CMS, ClassicalCMS, BayesianCMS
from cms.conformal import ConformalCMS
from cms.bootstrap import BootstrapCMS
from cms.diagnostics import evaluate_marginal, evaluate_conditional

# Parse input arguments
print ('Number of arguments:', len(sys.argv), 'arguments.')
print ('Argument List:', str(sys.argv))
model_num = 1
if len(sys.argv) != 14:
    print("Error: incorrect number of parameters.")
    quit()

#########################
# Experiment parameters #
#########################

# Default
sketch_name = "cms-cu"
stream_name = "dp-1"
d = 5
w = 20
n = 2000
method = "conformal"
method_unique = False
n_bins = 1
n_track = 1000
seed = 2021
posterior = "betamix"
confidence = 95/100
two_sided = False

# Input
sketch_name = sys.argv[1]
stream_name = sys.argv[2]
d = int(sys.argv[3])
w = int(sys.argv[4])
n = int(sys.argv[5])
method = sys.argv[6]
method_unique = sys.argv[7]
if method_unique != "NA":
    method_unique = bool(int(method_unique))
n_bins = sys.argv[8]
n_track = sys.argv[9]
seed = int(sys.argv[10])
posterior = sys.argv[11]
confidence_str = sys.argv[12]
two_sided = bool(int(sys.argv[13]))

confidence = float(confidence_str)/100

if n_bins != "NA":
    n_bins = int(n_bins)
if n_track != "NA":
    n_track = int(n_track)

# Stream parameters
n_docs = 2000

n_test=1000 # NOTE: originally 10000

# Output file
outfile_prefix = "exp2_" + sketch_name + "_" + stream_name + "_d" + str(d) + "_w" + str(w) + "_n" + str(n) + "_s" + str(seed) + "_" + posterior

##################
# Run experiment #
##################

np.random.seed(seed)

# # Define word data stream
# if stream_name == "words":
#     exchangeability = 1.0
#     filename = "data/words_nd" + str(n_docs) + "_exch" + str(exchangeability) + ".txt"
#     if os.path.exists(filename):
#         print("Loading stream from {:s}...".format(filename))
#         sys.stdout.flush()
#         stream = StreamFile(filename, seed=seed)
#     else:
#         print("Initializing stream from {:d} documents...".format(n_docs))
#         sys.stdout.flush()
#         stream = WordStream(n_docs=n_docs, exchangeability=exchangeability, filename=filename)
#     stream.set_seed(seed)
#     print("Loaded {:d} words.\n".format(len(stream.data)))
#     sys.stdout.flush()

# Define covid dna data stream
if stream_name == "covid":
    filename = "data/covid_dna.txt"
    if os.path.exists(filename):
        print("Loading stream from {:s}...".format(filename))
        sys.stdout.flush()
        stream = StreamFile(filename)
    else:
        print("Initializing stream...")
        sys.stdout.flush()
        stream = DNAStream(k_size = 16, seq_max=1000, filename_out=filename)
    stream.set_seed(seed)
    print("Loaded {:d} words.\n".format(len(stream.data)))
    sys.stdout.flush()

# Define word data stream
elif stream_name == "words":
    filename = "data/words.txt"
    if os.path.exists(filename):
        print("Loading stream from {:s}...".format(filename))
        sys.stdout.flush()
        stream = StreamFile(filename)
    else:
        print("Initializing stream...")
        sys.stdout.flush()
        stream = WordStream(n_docs=100, n_grams=2, filename_out=filename)
    stream.set_seed(seed)
    print("Loaded {:d} words.\n".format(len(stream.data)))
    sys.stdout.flush()

elif stream_name.startswith("zipf-"):
    # Dirichlet process
    alpha = float(stream_name.split("-")[1])
    sigma = None
    tau = None
    stream = Zipf(alpha, seed=seed)

elif stream_name.startswith("dp-"):
    # Dirichlet process
    alpha = float(stream_name.split("-")[1])
    sigma = None
    tau = None
    stream = DP(alpha, seed=seed)

elif stream_name.startswith("sp-"):
    # Dirichlet process
    sigma = float(stream_name.split("-")[1])
    alpha = None
    tau = None
    stream = SP(sigma, seed=seed)

elif stream_name.startswith("pyp-"):
    # Dirichlet process
    alpha = float(stream_name.split("-")[1])
    sigma = float(stream_name.split("-")[2])
    tau = None
    stream = PYP(alpha, sigma, seed=seed)

elif stream_name.startswith("nggp-"):
    # NGGP
    alpha = float(stream_name.split("-")[1])
    sigma = float(stream_name.split("-")[2])
    tau = float(stream_name.split("-")[3])
    stream = NGGP(alpha, sigma, tau, seed=seed)

else:
    print("Unknown data stream name: {:s}.\n".format(stream_name))
    sys.stdout.flush()
    exit
    
# Initialize method
if sketch_name == "cms":
    cms = CMS(d, w, seed=seed, conservative=False)
elif sketch_name == "cms-cu":
    cms = CMS(d, w, seed=seed, conservative=True)
else:
    print("Error: unrecognized sketch name")
    quit()

if method == "conformal-bayesian-dp":
    worker = ConformalCMS(stream, cms,
                          n_track = n_track,
                          unique = method_unique,
                          n_bins = n_bins,
                          scorer_type = "Bayesian-DP",
                          two_sided=two_sided)
    method_name = method + "_unique" + str(int(method_unique)) + "_bins" + str(n_bins) + "_track" + str(n_track)

elif method == "conformal-bootstrap":
    worker = ConformalCMS(stream, cms,
                          n_track = n_track,
                          unique = method_unique,
                          n_bins = n_bins,
                          scorer_type = "Bootstrap",
                          two_sided=two_sided)
    method_name = method + "_unique" + str(int(method_unique)) + "_bins" + str(n_bins) + "_track" + str(n_track)

elif method == "conformal-constant":
    worker = ConformalCMS(stream, cms,
                          n_track = n_track,
                          unique = method_unique,
                          n_bins = n_bins,
                          scorer_type = "Constant", 
                          two_sided=two_sided)
    method_name = method + "_unique" + str(int(method_unique)) + "_bins" + str(n_bins) + "_track" + str(n_track)

elif method == "conformal-proportional":
    worker = ConformalCMS(stream, cms,
                          n_track = n_track,
                          unique = method_unique,
                          n_bins = n_bins,
                          scorer_type = "Proportional",
                          two_sided=two_sided)
    method_name = method + "_unique" + str(int(method_unique)) + "_bins" + str(n_bins) + "_track" + str(n_track)

elif method == "conformal-adaptive":
    worker = ConformalCMS(stream, cms,
                          n_track = n_track,
                          prop_train = 0.5,
                          unique = method_unique,
                          n_bins = n_bins,
                          scorer_type = "Adaptive", 
                          two_sided=two_sided)
    method_name = method + "_unique" + str(int(method_unique)) + "_bins" + str(n_bins) + "_track" + str(n_track)

elif method == "bootstrap":
    worker = BootstrapCMS(stream, cms, two_sided=two_sided)
    method_name = method

elif method == "classical":
    worker = ClassicalCMS(stream, cms)
    method_name = method

elif method == "bayesian-oracle":
    if stream_name.startswith("dp-"):
        model = "DP"
    elif stream_name.startswith("sp-"):
        model = "NGGP"
        alpha = 1
        tau = 0
    elif stream_name.startswith("nggp-"):
        model = "NGGP"
    else:
        print("Error! Unknown model.")
        pdb.set_trace()
    worker = BayesianCMS(stream, cms, model=model, alpha=alpha, sigma=sigma, tau=tau, posterior=posterior, two_sided=two_sided)
    method_name = method

elif method.startswith("bayesian-dp-oracle"):
    worker = BayesianCMS(stream, cms, model="DP", alpha=alpha, posterior=posterior, two_sided=two_sided)
    method_name = method

elif method.startswith("bayesian-dp-"):
    alpha_cms_str = method.split("-")[2]
    if alpha_cms_str == "?":
        alpha_cms = None
    else:
        alpha_cms = float(alpha_cms_str)
    worker = BayesianCMS(stream, cms, model="DP", alpha=alpha_cms, posterior=posterior, two_sided=two_sided)
    method_name = method

elif method.startswith("bayesian-sp-oracle"):
    worker = BayesianCMS(stream, cms, model="NGGP", alpha=1, sigma=sigma, tau=1, posterior=posterior, two_sided=two_sided)
    method_name = method

elif method.startswith("bayesian-sp-"):
    sigma_cms_str = method.split("-")[2]
    if sigma_cms_str == "?":
        sigma_cms = None
    else:
        sigma_cms = float(sigma_cms_str)
    worker = BayesianCMS(stream, cms, model="NGGP", alpha=1, sigma=sigma_cms, tau=1, posterior=posterior, two_sided=two_sided)
    method_name = method

elif method == "bayesian":
    worker = BayesianCMS(stream, cms, posterior=posterior, two_sided=two_sided)
    method_name = method

else:
    print("Error: unrecognized method name")
    quit()

# Run experiment
results = worker.run(n, n_test, confidence=confidence, seed=seed)

# Header
def add_header(df):
    df["sketch"] = sketch_name
    df["data"] = stream_name
    df["d"] = d
    df["w"] = w
    df["method"] = method
    df["method-unique"] = method_unique
    df["posterior"] = posterior
    df["n_bins"] = n_bins
    df["n_track"] = n_track
    df["n"] = n
    df["seed"] = seed
    df["confidence"] = confidence
    df["two_sided"] = two_sided
    return df

################
# Save results #
################
outfile = "results/" + stream_name + "/detailed/" + outfile_prefix + "_" + method_name + "_" + confidence_str + "_ts" + str(int(two_sided)) + ".txt"
add_header(results).to_csv(outfile, index=False)
print("\nDetailed results written to {:s}\n".format(outfile))
sys.stdout.flush()

###############
# Diagnostics #
###############

s1 = evaluate_marginal(results, include_seen=True)
s2 = evaluate_marginal(results, include_seen=False)
s1u = evaluate_marginal(results, include_seen=True, unique=True)
s2u = evaluate_marginal(results, include_seen=False, unique=True)
summary_marginal = pd.concat([s1, s2, s1u, s2u])
print("Marginal summary:")
print(summary_marginal)
print()
sys.stdout.flush()

outfile = "results/" + stream_name + "/marginal/" + outfile_prefix + "_" + method_name + "_" + confidence_str + "_ts" + str(int(two_sided)) + ".txt"
add_header(summary_marginal).to_csv(outfile, index=False)
print("\nMarginal summary written to {:s}\n".format(outfile))
sys.stdout.flush()

s1 = evaluate_conditional(results, nbins=5, include_seen=True)
s2 = evaluate_conditional(results, nbins=5, include_seen=False)
s1u = evaluate_conditional(results, nbins=5, include_seen=True, unique=True)
s2u = evaluate_conditional(results, nbins=5, include_seen=False, unique=True)
summary_conditional = pd.concat([s1, s2, s1u, s2u])
print("Conditional summary:")
print(summary_conditional)
print()
sys.stdout.flush()

outfile = "results/" + stream_name + "/conditional/" + outfile_prefix + "_" + method_name + "_" + confidence_str + "_ts" + str(int(two_sided)) + ".txt"
add_header(summary_conditional).to_csv(outfile, index=False)
print("\nConditional summary written to {:s}\n".format(outfile))
sys.stdout.flush()

#pdb.set_trace()
