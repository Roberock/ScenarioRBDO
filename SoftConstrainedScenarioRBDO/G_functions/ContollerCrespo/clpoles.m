function lsfval=clpoles(p,stru)

lsfval=[];
limit=feval(stru.thisreq.limit,stru.thisreq.t);
for i=1:size(p,1)
    lsfval(i,:)=dompole(p(i,:),stru.gains)-limit(1,:);
end
lsfval=lsfval';

function out=dompole(p,gains)
plant=build_plant(p);
compensator=build_comp(gains);
Bp1=plant.B(:,1);
Bp2=plant.B(:,2);
Acl=[plant.A-Bp1*compensator.D*plant.C,Bp1*compensator.C;-compensator.B*plant.C,compensator.A];
Bcl=[Bp2;zeros(size(compensator.A,1),1)];
Ccl=[plant.C,zeros(1,size(compensator.A,1))];
Dcl=0;
clsys=ss(Acl,Bcl,Ccl,Dcl);
%pole(clsys)
out=max(real(pole(clsys)));
