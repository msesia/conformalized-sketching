import numpy as np
import pandas as pd
from collections import defaultdict
import pdb

def sum_dict(x, y):
    z = x
    for k in y.keys():
        x[k] += y[k]
    return z

def dictToList(d):
    reg_list = [[x] * d[x] for x in d.keys()]
    flat_list = [item for sublist in reg_list for item in sublist]
    return flat_list

def listToDict(v):
    d = defaultdict(lambda: 0)
    for x in v:
        d[x] += 1
    return d

def sort_dict(d):
    return dict(sorted(d.items(), key=lambda x: x[1], reverse=True))
