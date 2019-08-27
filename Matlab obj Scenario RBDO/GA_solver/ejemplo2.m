
clear all
close all
yy=-500:10:500;
xx=-500:10:500;
for i=1:size(yy,2)
   for j=1:size(xx,2)
      z(i,j)=objfun_schwefel([xx(j),yy(i)]);
      x(i,j)=xx(j);
      y(i,j)=yy(i);
   end
end
surfc(x,y,z,'LineStyle','none','FaceColor','interp')


