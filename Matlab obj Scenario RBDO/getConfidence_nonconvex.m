
function epsilon=getConfidence_nonconvex(k,beta,N)

% this function calculate the reliability epsilon for a support subsamble
% of dimension k
%
 
% INPUT
% N       = number of scenarios (i.e. samples, multi-dimensional data points);
% d       = number of optimization variables;
% k       = number of support constraints
% beta   = confdence parameter;

% OUTPUT
% epsilon = reliability parameter (possibly small);
%% Use example

%%
if k == N
    epsilon=1;
elseif k < N
    epsilon=1-(beta./(N.*nchoosek(N,k))).^(1/(N-k));
else
    error('number of support constraints k > number of scenarios N') 
end 
end