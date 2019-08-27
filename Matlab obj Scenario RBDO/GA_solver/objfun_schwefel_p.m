function valor=objfun_schwefel_p(x,param)
%
% Schwefel function
% 
% Function : f(X1,X2)=-x1*sin(sqrt(abs(x1)))-x2*sin(sqrt(abs(x2)))
% Search space constraints : -500<=xi<=500
%
% Solution: min(f)=f6(420.9687,420.9687)=-837.9658
%
% param: useless parameter only to show how to pass parameters
% using parameters for wasting time ;)

for i=1:param.n
    valor=param.p*param.p;
end
valor=-x(1).*sin(abs(x(1)).^0.5)-x(2).*sin(abs(x(2)).^0.5);


