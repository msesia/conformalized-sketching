import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import pdb

def evaluate_marginal(results, include_seen=False, unique=False):
    if unique:
        results = results.groupby("x").agg(lambda x: x.iloc[0]).reset_index()

    def nunique(x):
        return len(unique(x))

    if not include_seen:
        if 'seen' in results.columns:
            results = results[results["seen"]==False].copy()

    results["unique"] = unique
    results["coverage"] = (results["lower"] <= results["count"]).astype(float)
    results["length"] = results["upper"] - results["lower"] + 1
    results["gap"] = results["count"] - results["lower"]
    results["error-mean"] = np.sqrt(np.power(results["count"] - results["mean"],2))
    results["error-median"] = np.sqrt(np.power(results["count"] - results["median"],2))
    results["error-map"] = np.sqrt(np.power(results["count"] - results["mode"],2))
    results["observations"] = 1
    results["distinct"] = results["x"]
    summary = results.groupby(["method","unique"]).agg({'distinct':'nunique', 'observations': 'sum', 'coverage' : 'mean', 'length' : 'mean', 'gap':'mean',
                                                        'error-mean' : 'mean', 'error-median' : 'mean', 'error-map' : 'mean'})
    summary["include-seen"] = include_seen
    return summary.reset_index()

def evaluate_conditional(results, nbins=5, include_seen=False, unique=False):
    if unique:
        results = results.groupby("x").agg(lambda x: x.iloc[0]).reset_index()

    if not include_seen:
        if 'seen' in results.columns:
            results = results[results["seen"]==False].copy()

    num_unique = len(set(results["count"]))
    if nbins > num_unique:
        nbins = num_unique

    results["count-bin"], bin_endpoints = pd.qcut(np.array(results["count"]), nbins, retbins=True, labels=False, precision=0, duplicates="drop")
    print("Bin endpoints")
    print(bin_endpoints)
    bin_labels = [str(int(bin_endpoints[i]+1))+"-"+str(int(bin_endpoints[i+1])) for i in range(len(bin_endpoints)-1)]
    results["count-range"] = np.array(bin_labels)[np.array(list(results["count-bin"])).astype(int)]
    print(bin_labels)
    results["unique"] = unique
    results["coverage"] = (results["lower"] <= results["count"]).astype(float)
    results["length"] = results["upper"] - results["lower"] + 1
    results["gap"] = results["count"] - results["lower"]
    results["observations"] = 1
    results["distinct"] = results["x"]
    summary = results.groupby(["method", "unique", "count-bin", "count-range"]).agg({'coverage' : 'mean',
                                                                                     'length' : 'mean',
                                                                                     'gap':'mean',
                                                                                     'observations' : 'sum',
                                                                                     'distinct':'nunique'})
    summary["include-seen"] = include_seen
    return summary.reset_index()
