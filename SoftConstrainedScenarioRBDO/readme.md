#  Codes to perform Soft-Constrained Scenario RBDO: Complexity-based Modulation of Failure Probability Bounds #


run MAIN_RBDO_withSoftScenarioConstraints_Extended.m


%%%  %%% Program-1 (soft-constrained joint reliability requirements)

% min_{d\in \Theta , \zeta^{(i)}>0} \ lbrace J(d) +\rho \sum\limits_{i=1}^{N}  \zeta^{(i)}
% Such that: w(d,\delta{(i)} \leq \zeta^{(i)} \rbrace

% where
% \delta are the available scenarios (samples of the uncertain factors)
%  d\in \Thetaa is a design vector (e.g. fitting coefficients, tunable parameters etc) in a convex design set \Theta
% J(d) is a convex cost function
% $\rho$ is a parameter weighting the cost of violating constraints and
% w(d,\delta)=\max\limits_{j\in\{1,..,n_g \}| g_j(d,\delta) % w is a convex worst-case reliability performance function,
% n_g  is the number of individual reliability requirements defined by the performance functions g_j j=1,...,n_g

%%% %%%  Program-2 (soft-constrained individual reliability requirements)

% min_{d\in \Theta , \zeta_j^{(i)}>0} \ lbrace J(d) +\sum\limits_{j=1}^{n_g} \rho_j \sum\limits_{i=1}^{N}  \zeta_j^{(i)}
% Such that: g_j(d,\delta) \leq \zeta_j^{(i)} i=1,...,N,~j=1,..,n_g\rbrace
% where $\rho_j$ are parameters weighting the cost of violation on the reliability requirement g_j

% for this probelm the support scenarios (complexity0 will be equal to
% S_N^*= the number of active constriaints + the number of violating constraints