function g=G_roofTrust(d,delta)
% minimize  J=@(d) 20224*d(1)+364*d(2);%=20224*As+364*Ac
% s.t. LBd=[ 0.0006,0.018 ];UBd=[0.0012,0.063]; % 0.0006 ? As ? 0.0012, 0.018 ? Ac ? 0.063
% g=-(0.03-(q*l.^2).*(3.81./((Ac+ErrorAc).*Ec)+1.13./((As+ErrorAs).*Es)))
% and g>0 means failure
dn=[0.001, 0.042]; % nominal desing
As=d(1);
Ac=d(2);
q=delta(:,1);   %  mu=20000 ; std=1400; %normal random Load
l=delta(:,2);   %  mu=12 ; std=0.12; %normal random Length factors for truss members
Es=delta(:,3);  %  mu=10^11 ; std=6*10^9; %Elastic modulus of steel
Ec=delta(:,4);  %  mu=2*10^10 ; std=1.2*10^9; %Elastic modulus of concrete
ErrorAs=delta(:,5);%  mu=0 ; std=5.9853*10^(-5); %error on the cross-sectional area of steel bars
ErrorAc=delta(:,6);%  mu=0 ; std= 0.0048; %error on the cross-sectional area of concrete bars
MaxLoad=0.03;
%MaxLoad=0.02;
g= ((q.*l.^2)/2).*(3.81./((Ac-ErrorAc).*Ec)+1.13./((As-ErrorAs).*Es))-MaxLoad;
end