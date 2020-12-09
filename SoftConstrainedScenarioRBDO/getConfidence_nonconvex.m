
function epsilon=getConfidence_nonconvex(k,N,beta) 
% this function calculate the reliability epsilon for a support subsamble
% of dimension k 
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
    % this is numerically  for high N and k
 %  epsilon=1-(beta./(N.*nchoosek(N,k))).^(1/(N-k)); 
     %% This expansion is slower but guarantees stability (and numerical accuracy)
epsilon=(beta*k/N^2)^(1/(N-k));
for l=1:k-1
   Temp=(1/((N-l)/(k-l)))^(1/(N-k));
   epsilon=epsilon* Temp;
end
epsilon=1-epsilon; 
else
    error('number of support constraints k > number of scenarios N') 
end 
end