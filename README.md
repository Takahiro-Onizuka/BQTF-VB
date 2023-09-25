# BQTF-VB

This repository provides R code implementing calibrated variational Bayes approximation for Bayesian quantile trend filtering, as proposed by the following paper.

Onizuka, T., Hashimoto. S. and Sugasawa, S. (2022), Fast and Locally Adaptive Bayesian Quantile Smoothing using Calibrated Variational Approximations. *arXiv:2211.04666*.

The repository includes the following files.

* ```BQTF-function.R```: Implementation of Gibbs sampling and variational Bayes approximation with/without calibration for Bayesian quantile trend filtering.

* ```BQTF-example.R```: One-shot example of fitting Bayesian quantile trend filtering under horseshoe prior.

* ```true-function.R```: The true data-generating functions such as piecewise constant and varying smoothness.
