function g=g_LinearCompetitive(delta,d)

%% see also Rocchetta at al. https://doi.org/10.1016/j.ress.2019.106755
% input
% d = desing vector
% delta = matrix of uncertain parameters  [Number of samples x Number of uncertain factors]
% output
% g = realizations of the reliability performance function 
 
%% G1: 
g1=delta(:,2)./d(1) + delta(:,1)./d(2) - d(3);
%% G2:  
g2= d(1)*(delta(:,1)) - delta(:,2)./d(2) - d(3);
%% combine
g=[g1 g2]; %g>0 is failure 
end