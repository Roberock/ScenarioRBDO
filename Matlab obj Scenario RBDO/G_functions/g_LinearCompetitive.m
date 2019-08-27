function g=g_LinearCompetitive(delta,d)
%% G1: 
g1=-1/d(1)*(delta(:,2))-1/d(2)*(delta(:,1))+d(3);
%% G2: 
%g2=d(4)-delta(:,2)-d(5).*delta(:,1).^2+d(6).*delta(:,1).^3; 
g2=+d(1)*(delta(:,1))-1/d(2)*(delta(:,2))-d(3);
%%
g=[-g1 g2]; %g>0 is failure
%% G3: Highly non-linear LS
 
end