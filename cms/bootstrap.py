import numpy as np
import sys
import pandas as pd
from scipy.stats.mstats import mquantiles
from tqdm import tqdm
import pdb

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
    if (a_set & b_set):
        return True
    else:
        return False

class BootstrapCMS:
    def __init__(self, stream, cms):
        self.stream = stream
        self.cms = cms

    def run(self, n, n_test, confidence=0.9, seed=2021):
        print("Running bootstrap method with n = {:d}...".format(n))
        sys.stdout.flush()

        # Process stream
        for i in tqdm(range(n)):
            x = self.stream.sample()
            self.cms.update_count(x)

        # CMS parameters
        r,w = self.cms.count.shape

        def estimate_noise_dist(x):
            n_mc = 10000
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

        # Evaluate
        print("Evaluating on test data....")
        sys.stdout.flush()
        np.random.seed(seed)
        results = []
        noise = None
        for i in tqdm(range(n_test)):
            x = self.stream.sample()
            y = self.cms.true_count[x]
            upper = self.cms.estimate_count(x)
            if noise is None:
                noise = estimate_noise_dist(x)
            delta = mquantiles(noise, confidence)[0]
            lower = np.maximum(0, upper - delta)

            # Estimation
            delta_median = mquantiles(noise, 0.5)[0]
            est_median = np.maximum(0, upper - delta_median)

            results.append({'method':'Bootstrap',
                            'x':x, 'count':y, 'upper': upper, 'lower':lower,
                            'mean':est_median, 'median':est_median, 'mode':est_median,
                            'seen':0})

        results = pd.DataFrame(results)
        results = results.sort_values(by=['count'], ascending=False)
        return results
