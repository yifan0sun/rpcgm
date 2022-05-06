# RP-CGM

## Code for use
Reweighted Penalized Conditional Gradient Method

Use example.ipynb to create an example problem and run method + screening.


 - problem.py 
    - initiates a task, exposes gradient and function evaluations and L-smoothness factor w.r.t. 1-norm
    
 - gauge.py 
    - contains gauges, exposes LMO, gauge evaluations, dual gauge evaluations, screening. 
    
 - scalarfn.py
   - contains gamma and phi constructions. Exposes  max and min slopes, function and gradient evaluations
   
 - method.py
   - encodes P-CGM and RP-CGM. 

## Code that reproduces figures in paper

 - small_experiments: small sensing experiment, using Gaussian random sensing matrix
 - groupnorm/nucnorm/TVnorm: shows method behavior for interesting norms
 - Dorothea: sparse logistic regression task. dataset (https://archive.ics.uci.edu/ml/datasets/dorothea)

## Paper:
Yifan Sun and Francis Bach. "Screening for a Reweighted Penalized Conditional Gradient Method." (Under review.) arXiv preprint arXiv:2107.01106 (2021).
Link: https://arxiv.org/abs/2107.01106
