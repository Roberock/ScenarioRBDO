function delta=DGM_Betasquared_rebuttal(N,Ndelta)
Nbetas=ceil(Ndelta/2);
a=0.1*ones(1,Nbetas);
b=2*ones(1,Nbetas);
delta=zeros(N,Ndelta);

for i=1:N 
    delta(i,1:Nbetas)= betarnd(a,b)+ betarnd(b,a).^2-rand();
    delta(i,Nbetas+1:end) = normrnd(zeros(1,Ndelta-Nbetas),ones(1,Ndelta-Nbetas)).^3;
end
end