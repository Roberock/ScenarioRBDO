function g=evaluate_g2(p,d,do_linear,tvec,show)

plant=build_plant(p);
gains.num=d(1:4);
gains.den=d(5:end);
sscom=build_comp(gains);
Bp1=plant.B(:,1);
Bp2=plant.B(:,2);
com_poles=eigs(sscom.A);
Acl=[plant.A-Bp1*sscom.D*plant.C,Bp1*sscom.C;-sscom.B*plant.C,sscom.A];
Bcl=[Bp2;zeros(size(sscom.A,1),1)];
Ccl=[plant.C,zeros(1,size(sscom.A,1))];
Dcl=0;
clsys=ss(Acl,Bcl,Ccl,Dcl);
numclpoles=length(pole(clsys));
poles=[pole(clsys);com_poles]; %poles of closed loop and of controller
g(1)=max(real(poles));  %hurwitz stability
%[poles,gm,pm]=margins(gains,p,0);
if do_linear
    [y,t]=impulse(clsys,tvec);
       if any(isnan(y)) || any(isinf(y))
         error('NaN or Inf in the y impulse')
    end
    u=lsim(ss(sscom.A,sscom.B,sscom.C,sscom.D),y,t);    
else
    [ylin,tlin,xlin]=impulse(clsys,tvec);
    way=2;
    dist=0*tvec;
    if way==1
        dist(2)=1;
        ic=zeros(9,1);
    else
    	ic=xlin(2,:); %equate ic from linsim after impulse
    end
    [t,state]=ode45(@(t,y)difeq(t,y,p,sscom,dist,tvec),tvec,ic);
    y=state(:,2);
    u=-state(:,5);
    %plot(t,y);hold on;plot(tlin,ylin,'r--');input('u')
end
pos=t>15;
maxy=0.1;
g(2)=max(abs(y(pos)))-maxy;
maxu=0.5;
g(3)=max(abs(u))-maxu;

if show
    set(figure,'Position',[100 400 1200 600]);
    trans=1;
    if trans
        rpoles=exp(real(poles))-1;
    else
        rpoles=real(poles);
    end
    tem=max(abs(imag(poles)));
    subplot(131);plot(rpoles(1:numclpoles),imag(poles(1:numclpoles)),'*');hold on;
    plot(rpoles(numclpoles+1:end),imag(poles(numclpoles+1:end)),'ks');
    plot([0,0],tem*[-1,1],'r--');axis([min(rpoles),max([0.1,max(rpoles),-min(rpoles)]),-tem,tem]);
    if trans==0
        xlabel('Re');
    else
        xlabel('e^{Re}-1');
    end
    ylabel('Im');set(gca,'FontName','Arial','FontSize',20);title(['g_1=',num2str(g(1))]);
    subplot(132);plot(t,y,'LineWidth',2);hold on;plot([15,25],maxy*[1,1],'r--');plot([15,25],-maxy*[1,1],'r--');axis([0,25,min([-0.2,min(y)]),1.1*max(y)]);
    xlabel('Time');ylabel('x_2');set(gca,'FontName','Arial','FontSize',20);title(['g_2=',num2str(g(2))]);
    val=max([maxu+0.1,max(abs(u))]);
    subplot(133);plot(t,u,'LineWidth',2);hold on;plot(t,u*0+maxu,'r--');plot(t,u*0-maxu,'r--');axis([0,tvec(end),-val,val])
    xlabel('Time');ylabel('u');set(gca,'FontName','Arial','FontSize',20);title(['g_3=',num2str(g(3))]);
end


