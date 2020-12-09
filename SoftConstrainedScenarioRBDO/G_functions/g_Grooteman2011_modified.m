function g=g_Grooteman2011_modified(delta,d)

%% see also Rocchetta at al. https://doi.org/10.1016/j.ress.2019.106755
% input
% d = desing vector
% delta = matrix of uncertain parameters  [Number of samples x Number of uncertain factors]
% output
% g = realizations of the reliability performance function 
 
%% G1: 
g1=-d(1)+ delta(:,1)+5*d(2)* delta(:,2) -2*d(3)*(delta(:,1)-delta(:,2)).^2;%LSF linear dependency with d 
%% G2:  
g2=-d(1)*(1-delta(:,2))+d(2).*delta(:,1).^2-d(3).*delta(:,1).^3; %LSF linear in d saddle point/concave with delta

g=[g1 g2]; %g>0 is failure
end