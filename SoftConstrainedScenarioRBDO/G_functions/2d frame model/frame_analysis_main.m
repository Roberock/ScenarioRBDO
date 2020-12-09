clear all
clc

tic

%% giving input

input_data = 'frame_input.xlsx';
ndata = xlsread(input_data, 'Nodal Information');
mdata = xlsread(input_data, 'Member Information');

%% connectivity

no = length(ndata(:,1));
nm = length(mdata(:,1));

n1 = mdata(:,2);
n2 = mdata(:,3);

x1 = ndata(n1,2);
y1 = ndata(n1,3);
x2 = ndata(n2,2);
y2 = ndata(n2,3);

%% material and geometric properties

l = sqrt((x2-x1).^2 + (y2-y1).^2);
A = mdata(:,4);
I = mdata(:,5);
E = mdata(:,6);
G = mdata(:,7);
k = mdata(:,8);
beta = (E.*I./(l.^3))./(k.*G.*A./l);
a = mdata(:,9);
b = mdata(:,10);

%% loading information

py = mdata(:,11); 
ay = mdata(:,12);
q1 = mdata(:,13); 
q2 = mdata(:,14); 
a1 = mdata(:,15); 
a2 = mdata(:,16); 
m = mdata(:,17); 
am = mdata(:,18);

%% initialisation

dof= zeros(6,1,nm);
T = zeros(6,6,nm);
K_e_local= zeros(6,6,nm);
K_e_global= zeros(6,6,nm);
Fj_e_local= zeros(6,1,nm);
Fj_e_global= zeros(6,1,nm);
K_assembled = zeros(3*no,3*no);
Fj_assembled = zeros(3*no,1);
D_member_local= zeros(6,1,nm);
F_member= zeros(6,1,nm);

%% obtaining transformation & stiffness matrices, joint load vectors and their assembly

for i = 1:nm
    
    T(:,:,i) = transformation_matrix(l(i),x1(i),y1(i),x2(i),y2(i));
    
    K_e_local(:,:,i) = stiffness_matrix(l(i),A(i),I(i),E(i),beta(i));
    K_e_global(:,:,i) = T(:,:,i)'*K_e_local(:,:,i)*T(:,:,i);
    
    Fj_e_local(:,:,i) = equivalent_joint_load(l(i),py(i),ay(i),q1(i),q2(i),a1(i),a2(i),m(i),am(i));
    Fj_e_global(:,:,i) = T(:,:,i)'*Fj_e_local(:,:,i);
    
    dof(:,:,i) = [3*n1(i)-2; 3*n1(i)-1; 3*n1(i); 3*n2(i)-2; 3*n2(i)-1; 3*n2(i)];
    K_assembled(dof(:,:,i),dof(:,:,i)) = K_assembled(dof(:,:,i),dof(:,:,i)) + K_e_global(:,:,i);
    Fj_assembled(dof(:,:,i),1) = Fj_assembled(dof(:,:,i),1) + Fj_e_global(:,:,i);
    
end

%% considering support settlement & nodal load

for i = 1:no
    
    support_settlement(3*i-2:3*i,1) = ndata(i,4:6)';
    F_nodal(3*i-2:3*i,1) = ndata(i,7:9)';
    restraints(3*i-2:3*i,1) = ndata(i,10:12)';
    
end

F_total = F_nodal - Fj_assembled;

%% identifying restraining degrees of freedom

total_dof = (1:3*no)';
restrained_dof = find(restraints == 1);
unrestrained_dof = setdiff(total_dof,restrained_dof);

%% getting active & restrained part of stiffness matrix & force vector

Kaa = K_assembled(unrestrained_dof,unrestrained_dof);
Kra = K_assembled(restrained_dof,unrestrained_dof);
Kar = K_assembled(unrestrained_dof,restrained_dof);
Krr = K_assembled(restrained_dof,restrained_dof);

F_active = F_nodal(unrestrained_dof);
F_fixed_active = Fj_assembled(unrestrained_dof);
F_fixed_restrained = Fj_assembled(restrained_dof);

D_restrained = support_settlement(restrained_dof);

%% solving displacement

D_active = Kaa\((F_active-F_fixed_active) - Kar*D_restrained);

F_restrained =  F_fixed_restrained + Kra*D_active + Krr*D_restrained;

D = support_settlement;
D(unrestrained_dof) = D_active;

%% obtaining member forces

for i = 1:nm
    
    D_member_local(:,:,i) = T(:,:,i)*D(dof(:,:,i));
    F_member(:,:,i) = Fj_e_local(:,:,i) + K_e_local(:,:,i)*D_member_local(:,:,i);
    
end

%%

toc