# RP-CGM
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
