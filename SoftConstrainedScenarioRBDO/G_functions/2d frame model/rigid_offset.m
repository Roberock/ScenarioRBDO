function [T_rigid_offset] = rigid_offset(a,b,l)

k1 = [1 0 0; 0 1 a*l; 0 0 1]; 
k3 = [1 0 0; 0 1 -b*l; 0 0 1]; 
k2 = zeros(3);

T_rigid_offset = [k1 k2; k2 k3];

end