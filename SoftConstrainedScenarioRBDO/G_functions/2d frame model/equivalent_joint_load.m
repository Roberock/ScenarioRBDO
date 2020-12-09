function [Fj] = equivalent_joint_load(l,py,ay,q1,q2,a,b,m,am)

if py ~= 0
    f31 = py*ay*(1-ay/l)^2 ;
    f61 = -py*((ay/l)^2)*(l-ay);
    f21 = py*(1 - ay/l) + (f31 + f61)/l;
    f51 = py*ay/l - (f31 + f61)/l;
    Fj1 = [0;f21;f31;0;f51;f61];
else
    Fj1 = zeros(6,1);
end

if q1 ~= 0 && q2 ~= 0
    c = @(x) (q1 + (q2 - q1)*(x-a)./b).*x.*(1-x./l).^2;
    d = @(x) (q1 + (q2 - q1)*(x-a)./b).*(l-x).*(x./l).^2;
    e = @(x) (q1 + (q2 - q1)*(x-a)./b).*(1-x./l);
    f = @(x) (q1 + (q2 - q1)*(x-a)./b).*x./l;
    f32 = integral(c, a,b+a) ;
    f62 = -integral(d, a,b+a);
    f22 = integral(e, a,b+a) + (f32 + f62)/l;
    f52 = integral(f, a,b+a) - (f32 + f62)/l;
    Fj2 = [0;f22;f32;0;f52;f62];
else
    Fj2 = zeros(6,1);
end

if m ~= 0
    f33 = -m*(1-am/l)*(1-3*am/l);
    f63 = m*am*(2*l-3*am)/l^2;
    f23 = (f33 + f63)/l;
    f53 = -(f33 + f63)/l;
    Fj3 = [0;f23;f33;0;f53;f63];
else
    Fj3 = zeros(6,1);
end

Fj = (Fj1 + Fj2 + Fj3);

end