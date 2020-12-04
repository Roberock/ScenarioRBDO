# Given-data-Reliability-based-Desing-Optimization-via-Scenario-Theory
## LIST OF METHODS FOR THE CLASS ScenarioRBDO %%%%%%

This class introduces a set of methods and proprieties to perform reliability-based-design-optimization by Scenario theory. 
Scenario optimization makes direct use of the available data (the uncertain parameters delta) 
thereby eliminating the need for estimating the distribution of the uncertain parameters.

Furthermore, scenario theory enables rigorously quantifying the probability of the resulting design satisfying the reliability requirements
imposed upon it regarding future, unseen data. (see Robustness methods) 


%  For further reading:
%  REFERENCES
% 'A Scenario Optimization Approach to Reliability-Based Design'  R.Rocchetta, L.G. Crespo, S.P. Kenny......


%  Reliability Methods: 

% INPUT: a system desing (d)

%      Compute_Gfun(d): evaluates the performance function g=[g1,..,gNg]
%      Compute_FailureProbability(d): evaluates the overall Pf on the available scenarios
%      Compute_maxG(d):
%      Compute_W(d)
%      Compute_ReliabilityMetrics(design)

%  Scenario Optimization Methods:
 
%           Th program SP2 minimize Pf (to this end it does not induce constraints)
%   NLCon-  NonLinearConstraint(theta,alpha,Gexamined,Gdeltaidx) % evaluate
%           non linear constraints for program SP1 or SP3

%  SP-1)    Optimize_SP1(alpha,Theta0):  minimizes alpha percentile of w (using fmincon)
%  SP-2)    Optimize_SP2(): minimizes Pf given-data (using GA);
%  SP-3)    Optimize_SP3(alpha,Theta0,Gexamined,Gdeltaidx):minimizes the sum
%           of the alpha percentiles of each requitement g_j with j=1,..,Ng  (using fmincon)

% Robustness method
%   getEpsilon(k,beta)  gets the non-convex robustness given cardinality k  and confidence beta (scenario size N is within the object)
 
%   ScenarioConstraints_addMethod (for SP1 and SP3)
%   ScenarioConstraints_removeMethod (for SP1 and SP3)

%  Outlier Removal Method
%   RemoveConstraints re-optimizes removing a list of scenario from the initial data set

%  Plots
%    plot_deta_vs_G scatter: the scenarios in the uncertainty space (2-D) vs the performance function realizations
%    plot_detaIndex_vs_G sort and plot the gj and the scenarios indices w.r.t. one of the reliability requirement
%    plot_2D_SafeFailDomains_and_Scenarios % plot failure and safe regions
%    Plot_ScenarioConstraints : plot a list of scenario constraints
