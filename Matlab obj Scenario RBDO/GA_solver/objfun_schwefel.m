function value=objfun_schwefel(x)
%
% Schwefel function
% 
% Function : f(X1,X2)=-x1*sin(sqrt(abs(x1)))-x2*sin(sqrt(abs(x2)))
% Search space constraints : -500<=xi<=500
%
% Solution: min(f)=f6(420.9687,420.9687)=-837.9658
%

value=-x(1).*sin(abs(x(1)).^0.5)-x(2).*sin(abs(x(2)).^0.5);