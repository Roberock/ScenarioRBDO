# ScenarioRBDO
 
## Soft-Constrained Scenario RBDO: Complexity-based Modulation of Failure Probability Bounds


<p align="center">
  <img src="./figs/SoftConstrainedExample.png.png" width="738">
</p>


Codes to prun the analysis 
**run MAIN_RBDO_withSoftScenarioConstraints_Extended.m

### Program-1 (soft-constrained joint reliability requirements)
```
 % min_{d\in \Theta , \zeta^{(i)}>0} \ lbrace J(d) +\rho \sum\limits_{i=1}^{N}  \zeta^{(i)}
 Such that: w(d,\delta{(i)} \leq \zeta^{(i)} \rbrace

 where
 \delta are the available scenarios (samples of the uncertain factors)
 d\in \Theta is a design vector (e.g. fitting coefficients, tunable parameters etc) in a convex design set \Theta
 J(d) is a convex cost function
 $\rho$ is a parameter weighting the cost of violating constraints and
 w(d,\delta)=\max\limits_{j\in\{1,..,n_g \}| g_j(d,\delta) % w is a convex worst-case reliability performance function,
 n_g  is the number of individual reliability requirements defined by the performance functions g_j j=1,...,n_g
 ```
### Program-2 (soft-constrained individual reliability requirements)
```
% min_{d\in \Theta , \zeta_j^{(i)}>0} \ lbrace J(d) +\sum\limits_{j=1}^{n_g} \rho_j \sum\limits_{i=1}^{N}  \zeta_j^{(i)}
% Such that: g_j(d,\delta) \leq \zeta_j^{(i)} i=1,...,N,~j=1,..,n_g\rbrace
% where $\rho_j$ are parameters weighting the cost of violation on the reliability requirement g_j

% for this probelm the support scenarios (complexity is S_N*)
% S_N^*= the number of active constriaints + the number of violating constraints
```
### Reliability bounds via Scenario optimization theory
Pobability[Pf(d*)<\epsilon]>1-\beta
where \beta is a small confidence parameter selected by the analyst and
\epsilon is a reliabiity bounds provided by scenario theory
fixing a confidence \epsilon is a function f(N,\beta,sN^*)
where N is the number of samples in program 1 and 2 and
s_N^* start is the complexity of the solution

### References

**generalization bounds for convex scenario program with relaxed constraints 

[] Garatti, Simone & Campi, Marco. (2019). Risk and complexity in scenario optimization. Mathematical Programming. https://doi.org/10.1007/s10107-019-01446-4

**Scenario RBDO for non-convex problems

[] Roberto Rocchetta, Luis G. Crespo, Sean P. Kenny, A scenario optimization approach to reliability-based design, Reliability Engineering & System Safety, Volume 196, 2020, 106755, ISSN 0951-8320, https://doi.org/10.1016/j.ress.2019.106755.

 









### Description of the CLASS ScenarioRBDO  (scenario RBDO for non convex problems)

This class introduces a set of methods and proprieties to perform reliability-based-design-optimization by Scenario theory. 
Scenario optimization makes direct use of the available data (the uncertain parameters delta) 
thereby eliminating the need for estimating the distribution of the uncertain parameters.

Furthermore, scenario theory enables rigorously quantifying the probability of the resulting design satisfying the reliability requirements
imposed upon it regarding future, unseen data. (see Robustness methods) 

 


### Reliability Methods: 

% INPUT: a system desing (d)

        Compute_Gfun(d): evaluates the performance function g=[g1,..,gNg]
        Compute_FailureProbability(d): evaluates the overall Pf on the available scenarios
        Compute_maxG(d):
        Compute_W(d)
        Compute_ReliabilityMetrics(design)

###  Scenario Optimization Methods:
 
     Th program SP2 minimize Pf (to this end it does not induce constraints)
     NLCon-  NonLinearConstraint(theta,alpha,Gexamined,Gdeltaidx) % evaluate
     non linear constraints for program SP1 or SP3

    SP-1)    Optimize_SP1(alpha,Theta0):  minimizes alpha percentile of w (using fmincon)
    SP-2)    Optimize_SP2(): minimizes Pf given-data (using GA);
    SP-3)    Optimize_SP3(alpha,Theta0,Gexamined,Gdeltaidx):minimizes the sum
             of the alpha percentiles of each requitement g_j with j=1,..,Ng  (using fmincon)

###  Robustness method
     getEpsilon(k,beta)  gets the non-convex robustness given cardinality k  and confidence beta (scenario size N is within the object)
 
     ScenarioConstraints_addMethod (for SP1 and SP3)
     ScenarioConstraints_removeMethod (for SP1 and SP3)

###   Outlier Removal Method
     RemoveConstraints re-optimizes removing a list of scenario from the initial data set

 
###   Data visualization methods
      plot_deta_vs_G scatter: the scenarios in the uncertainty space (2-D) vs the performance function realizations
      plot_detaIndex_vs_G sort and plot the gj and the scenarios indices w.r.t. one of the reliability requirement
      plot_2D_SafeFailDomains_and_Scenarios % plot failure and safe regions
      Plot_ScenarioConstraints : plot a list of scenario constraints
