function g = Lateral_motion_Controller_statespace(A,B,W,P)
 

tvec=0:0.01:20;
K=W'*P;

com_poles=eigs(P);
 
clsys = ss(A,B,0,0);

numclpoles=length(pole(clsys));


numclpoles=length(pole(clsys));
poles=[pole(clsys);com_poles]; %poles of closed loop and of controller
g(1)=max(real(poles));  %hurwitz stability


[y,t]=impulse(clsys,tvec);
  u=lsim(clsys,y,t);    
  
pos=t>15;
maxy=0.1;
g(2)=max(abs(y(pos)))-maxy;
maxu=0.5;
g(3)=max(abs(u))-maxu;