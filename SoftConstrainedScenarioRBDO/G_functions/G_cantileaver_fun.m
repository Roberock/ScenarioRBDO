function [g] = G_cantileaver_fun(delta,d)
 
% d=(b,h,L) geometry of the beam
% P=(Px,Py,S,E,Error_length) meterial propriety and loads 
Max_deflection = 2.5; % the max deflection
b = d(1); 
h = d(2); 
Ldet = d(3); 
% uncertain factors
Px=delta(:,1); % random force  in x     : Normal (Mean,std)  (500, 100);
Py=delta(:,2); % random force  in y     : Normal (Mean,std)  (1000, 100);
S=delta(:,3);  % the max yield stress   : Normal (Mean,std)  (4.0e4, 2.0e3);
E=delta(:,4);  % Young’s Modulus        : Normal (Mean,std)  (29.0e6, 3.0e6);
errL=delta(:,5);  % Error on the length : Normal (Mean,cv)   (0, 0.01);
% N=10^5;
%  delta(:,[1,2])=mvnrnd([500 10^3],[100 90; 90 150],N); % assume Px and Py are correlated
%  delta(:,3)=normrnd(4.0e4,2.0e3,[N,1]); %  
%  delta(:,4)=normrnd(29.0e6, 3.0e6,[N,1]); % 
%  delta(:,5)=normrnd(0,0.01,[N,1]); % error on the length of the bar 

L=Ldet.*errL;
%Calculate g-performance
g(:,1) = 6.*L./(b.*h).*(Px./b+Py./h)- S; % non-convex in d(1), d(2) and delta 
g(:,2) = 4*L.^3./E.*((Px/(b.^3.*h)).^2+(Py/(b.*h.^3)).^2).^0.5 -Max_deflection; % non-convex in d(1), d(2) and delta 
end