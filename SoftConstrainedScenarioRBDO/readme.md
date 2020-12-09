# Soft-Constrained Scenario RBDO: Complexity-based Modulation of Failure Probability Bounds #

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
```
Pobability[Pf(d*)<\epsilon]>1-\beta
where \beta is a small confidence parameter selected by the analyst and
\epsilon is a reliabiity bounds provided by scenario theory
fixing a confidence \epsilon is a function f(N,\beta,sN^*)
where N is the number of samples in program 1 and 2 and
s_N^* start is the complexity of the solution

 ```
 
References
 ###  convex scenario program with relaxed constraints 
Bounds derived and applied in:
 -  [] Garatti, Simone & Campi, Marco. (2019). Risk and complexity in scenario optimization. 
        Mathematical Programming. https://doi.org/10.1007/s10107-019-01446-4
 
###  non-convex Scenario RBDO  

  - [] Roberto Rocchetta, Luis G. Crespo, Sean P. Kenny, A scenario optimization approach to reliability-based design,
     Reliability Engineering & System Safety, Volume 196, 2020, 106755, ISSN 0951-8320, https://doi.org/10.1016/j.ress.2019.106755.

