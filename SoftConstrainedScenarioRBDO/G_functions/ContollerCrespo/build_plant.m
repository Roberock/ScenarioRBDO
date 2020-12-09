function plant=build_plant(p)
%linearization about zero (not about stable oscillator). nonlinear spring does not show up
if length(p)==2
    plant.p=[p(1),1,p(2),0,1,0.001];
else
    plant.p(1)=p(1);    %m1  
    plant.p(2)=p(2);    %m2
    plant.p(3)=p(3);    %klinear
    %plant.p(4)=p(4);   %knonlin
    plant.p(5)=p(5);    %control effectivenes
    plant.p(6)=p(6);    %time delay
end
ratio=1/plant.p(6);
plant.A=[0,0,1,0,0;....
         0,0,0,1,0;....
         -plant.p(3)/plant.p(1),plant.p(3)/plant.p(1),0,0,plant.p(5)/plant.p(1);....
          plant.p(3)/plant.p(2),-plant.p(3)/plant.p(2),0,0,0;....
          0,0,0,0,-ratio];
Bp1=[0,0,0,0,ratio]';
Bp2=[0,0,0,1/plant.p(2),0]';
plant.B=[Bp1,Bp2];
plant.C=[0,1,0,0,0];
plant.D=[0,0];
plant.x0=Bp2;     %simulating impulse: zero input and IC=B