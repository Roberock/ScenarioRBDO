function g = g_4D_parallelSystem(delta,d)
%% Reliability of a parallel system
%see ''Structural reliability analysis of multiple limit state functions using multi-input multi-output support vector machine'' 
% Hong-Shuang Li  An-Long Zhao1 and Kong Fah Tee  

%% see also Rocchetta at al. https://doi.org/10.1016/j.ress.2019.106755
% input
% d = desing vector
% delta = matrix of uncertain parameters  [Number of samples x Number of uncertain factors]
% output
% g = realizations of the reliability performance function 

%% Rocchetta at al. https://doi.org/10.1016/j.ress.2019.106755
%% G1:
g1 = +d(1).*(delta(:,1)-delta(:,2)).^2  - (delta(:,1)+ delta(:,2))./sqrt(d(2))  +d(5);
%% G2:
g2 = +d(1).*(delta(:,1)-delta(:,2)).^2  + (delta(:,1)+ delta(:,2))./sqrt(d(2))  +d(5); 
%% G3:
g3 = -d(3)*delta(:,1)+delta(:,2)  - d(4)*sqrt(d(2)); 
%% G4:
g4 = +d(3)*delta(:,1)-delta(:,2)  - d(4)*sqrt(d(2));
%%
g=[ g1  g2  g3  g4]; % we have that g>0 implies failure 
end
