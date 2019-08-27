function g=g_4D_parallelSystem(delta,d)
% Reliability of a parallel system
%see ''Structural reliability analysis of multiple limit state functions using multi-input multi-output support vector machine'' Hong-Shuang Li  An-Long Zhao1 and Kong Fah Tee
%Nominal Desing
% d(1)=0.1
% d(2)=2
% d(3)=1
% d(4)=3.5
% d(5)=3
%% G1-4
g1=d(1).*(delta(:,1)-delta(:,2)).^2-( delta(:,1)+ delta(:,2))./sqrt(d(2))+d(5);
g2= d(1).*(delta(:,1)-delta(:,2)).^2+( delta(:,1)+ delta(:,2))./sqrt(d(2))+d(5);

g3=d(3)*delta(:,1)-delta(:,2)+d(4)*sqrt(d(2)); 
g4= -d(3)*delta(:,1)+delta(:,2)+d(4)*sqrt(d(2));
%%
g=[ g1  g2  g3  g4]; %g>0 is failure
g=-g;
%% G3: Highly non-linear LS

end
