import numpy as np
from sklearn.model_selection import train_test_split
from scipy.stats.mstats import mquantiles
import copy

from sklearn.preprocessing import SplineTransformer
from sklearn.linear_model import QuantileRegressor
from sklearn.isotonic import IsotonicRegression

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()


import pdb
import matplotlib.pyplot as plt

class IQR:
    """
    Isotonic quantile regression model
    """
    def __init__(self, quantile, quantiles, seed, log_transform=False):
        self.r = robjects.r
        #path = os.path.dirname(__file__)
        #filename = path + "/ratio_uniforms.R"
        #self.r['source'](filename)
        self.r('''library(isodistrreg) ''')
        self.quantile = np.round(quantile,5)
        self.quantiles = quantiles
        self.seed = seed
        self.log_transform = log_transform

    def fit(self, X, Y, X_calib=None, Y_calib=None, debug=False):
        self.r('''
        fit <- NULL
        fit_idr = function(X,Y) {
            X <- as.numeric(X)
            Y <- as.numeric(Y)
            fit <<- isodistrreg::idr(Y,data.frame(X))
            return(0)
        }
        ''')

        idx_sort = np.argsort(X.flatten())
        X = X[idx_sort]
        Y = Y[idx_sort]

        r_fit_idr = robjects.r['fit_idr']
        # Detect outliers
        #idx_inliers=[]
        idx_inliers = np.where(X < np.mean(X) + 3 * np.std(X))[0]
        if len(idx_inliers)>0:
            X_in = X[idx_inliers]
            Y_in = Y[idx_inliers]
        else:
            X_in = X
            Y_in = Y
        self.X_train = X_in
        r_fit_idr(X_in,Y_in)

        # Fit the spline model to interpolate
        self.qr = QR(quantiles=[self.quantile], log_transform=self.log_transform)
        self.qr.fit(X,Y)

    def predict(self,X,quantiles=None):
        self.r('''
        predict_idr = function(X,quantile,seed) {
            X <- as.numeric(X)
            quantile <- as.numeric(quantile)
            predictions <- predict(fit, data = data.frame(X), interpolation="linear", seed=seed)
            qhat <- as.vector(qpred(predictions, quantiles = quantile))
            return(qhat)
        }
        ''')
        X = np.array(X)

        r_predict_idr = robjects.r['predict_idr']
        if quantiles is None:
            pred = np.array(r_predict_idr(X, self.quantile, self.seed))
        else:
            pred = np.array([r_predict_idr(x, quantiles, self.seed) for x in X])

        # Interpolate
        if True:
            pred_qr = self.qr.predict(X)
            idx_above = np.where(X>np.max(self.X_train))[0]
            #idx_above = np.arange(len(X))
            if len(idx_above)>0:
                if(len(pred.shape)>1):
                    for j in range(pred.shape[1]):
                        pred[idx_above,j] = pred_qr[idx_above,0]
                else:
                    pred[idx_above] = pred_qr[idx_above]
            # idx_below = np.where(X<np.min(self.X_train))[0]
            # if len(idx_below)>0:
            #     pred[idx_below] = np.min(pred)

        if len(self.X_train)<100:
            pred = pred_qr[:,0]

        # Trim within reasonable domain
        if(len(pred.shape)>1):
            for j in range(pred.shape[1]):
                pred[:,j] = np.maximum(0, pred[:,j])
                pred[:,j] = np.minimum(X, pred[:,j])
                # Make sure there is no quantile crossing
                pred = -np.sort(-pred,1)
        else:
                pred = np.maximum(0, pred)
                pred = np.minimum(X, pred)


        return pred


class QR:
    """
    Quantile regression model
    """
    def __init__(self, quantiles, degree=1, log_transform=False):
        self.degree = degree
        self.quantiles = np.round(quantiles,5)
        self.augmentation = CustomTransformer()
        self.n_quantiles = len(self.quantiles)
        self.log_transform = log_transform

    def fit(self, X, Y, X_calib=None, Y_calib=None, debug=False):
        X = X.flatten()

        # Transform the variable
        X_aug = self.augmentation.fit_transform(X)
        # Fit all the quantile regression models
        self.qr = [QuantileRegressor(quantile=q, solver='highs', alpha=0) for q in self.quantiles]
        for i in range(len(self.quantiles)):
            self.qr[i].fit(X_aug, Y)


    def predict(self, X, quantiles=None):
        X = np.array(X)
        n = len(X)
        X_aug = self.augmentation.transform(X)
        output = np.zeros((n,self.n_quantiles))
        for i in range(len(self.quantiles)):
            output[:,i] = self.qr[i].predict(X_aug)
            # Trim within reasonable domain
            output[:,i] = np.maximum(0, output[:,i])
            output[:,i] = np.minimum(X, output[:,i])
        output = output
        return output

class QRScores:
    """
    CQR for lower bound
    By definition, lower scores result in higher lower bounds (more liberal).
    Higher scores result in more conservative bounds.
    """
    def __init__(self, cms, confidence, seed, two_sided=False):
        self.cms = copy.deepcopy(cms)
        self.two_sided = two_sided
        alpha = 1.0 - confidence
        if self.two_sided:
            self.quantile_low = alpha
            self.quantile_upp = 1.0-alpha
            self.bbox_low = IQR(self.quantile_low,[self.quantile_low],seed)
            self.bbox_upp = IQR(self.quantile_upp,[self.quantile_upp],seed)
        else:
            self.quantile = 1.0-alpha
            self.bbox = IQR(self.quantile,[self.quantile],seed)

    def train(self, X, Y, X_calib=None, Y_calib=None):
        # Fit black-box model
        if self.two_sided:
            self.bbox_low.fit(X, Y, X_calib=X_calib, Y_calib=Y_calib)
            self.bbox_upp.fit(X, Y, X_calib=X_calib, Y_calib=Y_calib)
        else:
            self.bbox.fit(X, Y, X_calib=X_calib, Y_calib=Y_calib)

    def compute_score(self, x, y):
        "This score measures by how much we need to decrease the upper bound to obtain a valid lower bound"
        upper_max = self.cms.estimate_count(x)
        if self.two_sided:
            lower_raw = self.bbox_low.predict([upper_max])[0]
            upper_raw = self.bbox_upp.predict([upper_max])[0]
            lower = np.maximum(0, np.minimum(lower_raw, upper_raw))
            upper = np.minimum(upper_max, np.maximum(lower_raw, upper_raw))
            score = np.maximum(lower-y, y-upper)

        else:
            lower = self.bbox.predict([upper_max])[0]
            score = lower - y # Positive score means the lower bound is too high

        return score.astype(int), 0

    def predict_interval(self, x, tau, tau_u):
        upper_max = self.cms.estimate_count(x)
        if self.two_sided:
            lower_raw = self.bbox_low.predict([upper_max])[0]
            upper_raw = self.bbox_upp.predict([upper_max])[0]
            lower = np.maximum(0, np.minimum(lower_raw, upper_raw)-tau)
            upper = np.minimum(upper_max, np.maximum(lower_raw, upper_raw)+tau)

        else:
            lower = self.bbox.predict([upper_max])[0] - tau
            lower = np.maximum(0, lower).astype(int)
            upper = upper_max

        return lower, upper_max

class CustomTransformer:
    def __init__(self, n_knots=50):
        self.n_knots = n_knots

    def fit_transform(self, X):
        self.knots = mquantiles(np.unique(X), np.arange(1,1+self.n_knots)/(self.n_knots+1)).astype(int)
        X_t = self.transform(X)
        return X_t

    def transform(self, X):
        K = self.n_knots
        X_t = np.tile(X, (K+1,1)).T
        for k in range(K):
            X_t[:,k+1] = (X - self.knots[k]) * (X >= self.knots[k])
        return X_t
