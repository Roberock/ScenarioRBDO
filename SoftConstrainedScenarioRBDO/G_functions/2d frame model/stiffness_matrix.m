function [K_local] = stiffness_matrix(l,A,I,E,bt)

a = E*A/l;
b = 12*E*I/(l^3*(1+12*bt));
c = 6*E*I/(l^2*(1+12*bt));
d = 4*E*I*(1+3*bt)/(l*(1+12*bt));
e = 2*E*I*(1-6*bt)/(l*(1+12*bt));

K_local = [a 0 0 -a 0 0;0 b c 0 -b c;0 c d 0 -c e;-a 0 0 a 0 0;0 -b -c 0 b -c;0 c e 0 -c d];

end
