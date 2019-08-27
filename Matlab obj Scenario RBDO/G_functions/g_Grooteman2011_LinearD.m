function g=g_Grooteman2011_LinearD(delta,d)

% g1= d1-d2(x1+x2)+d3(x1-x2)^2
% g2= d1-x2-d2*x1^2+d3*x1^3
g1=d(1)- 3*d(2).*(delta(:,1)+delta(:,2))+d(3)*(delta(:,1)-delta(:,2)).^2;%LSF polinomial in d convex quadratic in delta
g2=2-delta(:,2)-d(2).*delta(:,1).^2-d(3).*delta(:,1).^3; %LSF linear in d saddle in delta

g=-[g1 g2]; %g>0 is failure
end