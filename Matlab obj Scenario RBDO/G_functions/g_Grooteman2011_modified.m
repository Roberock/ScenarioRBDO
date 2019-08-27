function g=g_Grooteman2011_modified(delta,d)
%% G1: 
g1=-d(1)+ delta(:,1)+5*d(2)* delta(:,2) -2*d(3)*(delta(:,1)-delta(:,2)).^2;%LSF linear dependency with d 
% g1=d(1)-(delta(:,1)+ delta(:,2))./sqrt(d(2))+d(3).^2*(delta(:,1)-delta(:,2)).^2;%LSF polinomial in d convex quadratic in delta
% g1=d(1)-(delta(:,1)+ delta(:,2))./sqrt(d(2))+exp(d(3))*(delta(:,1)-delta(:,2)).^2;%LSF arbitrary dependency with d, convex quadratic in delta
%% G2: 
%g2=d(4)-delta(:,2)-d(5).*delta(:,1).^2+d(6).*delta(:,1).^3; 
g2=-d(1)*(1-delta(:,2))+d(2).*delta(:,1).^2-d(3).*delta(:,1).^3; %LSF linear in d saddle point/concave with delta

g=[g1 g2]; %g>0 is failure
%% G3: Highly non-linear LS
%g2=d(4)-delta(:,2)-d(5).*delta(:,1).^2+d(6).*delta(:,1).^3; 
%g3=3-d(1)*delta(:,2)+(d(4)*delta(:,1)).^4; %LSF linear in d saddle point/concave with delta

%g=[g1 g2 g3]; %g>0 is failure
end