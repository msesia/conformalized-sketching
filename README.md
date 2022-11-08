# Conformalized Sketching

This repository contains a software implementation of the methods described in the papers:
 - "Conformalized Frequency Estimation from Sketched Data", by Matteo Sesia and Stefano Favaro (accepted in NeurIPS 2022)
 - "Conformal Frequency Estimation with Sketched Data under Relaxed Exchangeability", by Matteo Sesia, Stefano Favaro, and Edgar Dobriban (longer version of the above, preprint)


This repository also contains the scripts utilized to carry out the numerical experiments and data analyses described in the aforementioned papers.

## Contents:
 - cms/             A Python package (with an R backend for some operations) implementing the proposed conformalized sketching method, as well as the 3 benchmark approaches (Classical, Bayesian, and Bootstrap)
 - experiments/     Scripts for reproducing the numerical experiments and producing the figures shown in the paper. The bash script "submit_experiments.sh" can be used to launch all experiments.
 - experiments/data Real data sets used in Sections 4.2 and 4.3 of Sesia and Favaro (2022)

## Software requirements:
   - Python >= 3.8.10
   - R >= 4.0.3
   - A few additional standard Python and R libraries as detailed in the package.