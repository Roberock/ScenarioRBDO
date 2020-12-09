function g=G_fun_lateralmotioncontroller(d,delta)
%% Design of a controller for the lateral motion for an aircraft
% Refer to the article:
% [1] "Probabilistic and Randomized Tools for Control Design"
% DOI: 10.1201/b10384-76 Fabrizio Dabbene Fabrizio et al
 
% Uncertain aircraft parameters:    Lp   L?   Lr  g/V  Y?   N?   Np   N?   Nr   L?a  Y?r   N?r   N?a
% In this work:     delta (1)  (2)  (3)  (4)  (5)  (6)  (7) (8)  (9)  (10)  (11)  (12)  (13)      
% Mean values  -2.93 -4.75 0.78 0.086 -0.11 0.1 -0.042 2.601 -0.29 -3.91 0.035 -2.5335 0.31
% Perturbated of 10% arround the mean
 
%  Design variables W and P (controller matrices)
%  Pseq = [0.3075 -0.3164 -0.0973 -0.0188
%        -0.3164 0.5822 -0.0703 -0.0993
%        -0.0973 -0.0703 0.2277 0.2661
%        -0.0188 -0.0993 0.2661 0.7100]
 
% Wseq = [-0.0191 0.2733
% -0.0920 0.4325
% 0.0803 -0.3821
% 0.4496 -0.2032]

% In this script it correspond to a desing d_opt [Probabilistic Sequential Design]:
% d_opt=[0.3075 -0.3164 -0.0973 -0.0188 -0.3164 0.5822 -0.0703 -0.0993 ...
%       -0.0973 -0.0703 0.2277 0.2661 -0.0188 -0.0993 0.2661 0.7100 -0.0191 ...
%        0.2733 -0.0920 0.4325 0.0803 -0.3821 0.4496 -0.2032]
 
% In this script it correspond to a desing d_opt [Scenario design]:
%   d_opt=[0.1445 -0.0728 0.0035 0.0085 -0.0728 0.2192 -0.0078 -0.0174 ...
%         0.0035 -0.0078 0.1375 0.0604 0.0085 -0.0174 0.0604 0.1975 ...
%       0.0109 0.0908 7.2929 3.4846 0.0439 -0.0565 0.6087 -3.9182]

 
% State-space equation is:
% dx(t)/dt = A(delta)x(t) + B(delta)u(t)

% where
% x(t) = lateral motion vector x(1) is the bank angle, x(2) its derivative, x(3) is the sideslip angle, x(3) the yaw rate; 
% dx(t)/dt = derivate of x over time;
% u(t)= controller u(1) the rudder deflection and u(2) the aileron deflection.

% The control goal is to design a state feedback controller u(t) = Kx(t)
% such that real parts of eigenvalues of the closedloop system are smaller than  -? < 0
% this means that the controller robustly stabilizes the system guaranteeing a desired decay rate ? > 0

% A quadratic performance criterion for this condition is  [2]
% f(W,P,delta)= A(delta)*P + P*A'(delta) + B(delta)*W' + W*B'(delta) + 2*?*P 
% if G(W,P,delta) is negative semidefinite the condition for stability is sattisified
% see [2] S. Boyd, L. El Ghaoui, E. Feron, and V. Balakrishnan. Linear Matrix Inequalities in System and Control Theory. SIAM, Philadelphia, 1994

%% The optimization problem can be defined as follows:
% \min\limits_{P \succ 0, W} \{ Tr[W] : P \succeq \thetaI, \max\lambda(f(W,P,delta) \preceq 0 \}
%where P \in \mathbb{R}^{4\times 4} and W \in \mathbb{R}^{2\times 4}. \theta= 0.01 is a small positive number to ensure the positive de?niteness of P. The control gain K can be ?nally computed as K = Y P?1
%% START
P= [d(1:4); 
    d(5:8); 
    d(9:12);
    d(13:16)];
% where P must be positive definite....  check if all(eig(P)>0) ?
W= [d(17:18);
    d(19:20);
    d(21:22); 
    d(23:24)];

alpha=0.1; % [3] A Posteriori Probabilistic Bounds of Convex Scenario Programs with Validation Tests Chao Shang et al 7 Jan 2020 arXiv:1903.11734v2 
g=zeros(size(delta,1),1);
parfor i=1:size(delta,1)
    A=[0    ,  1     ,    0       ,    0
       0 ,delta(i,1)         ,delta(i,2) , delta(i,3)
       delta(i,4)            ,   0       , delta(i,5)                       ,   -1
       delta(i,6)*delta(i,4) , delta(i,7), delta(i,6)+delta(i,6)*delta(i,5), delta(i,9)-delta(i,8)];

    B=[ 0   ,  0
        0         ,  delta(i,10)
        delta(i,11) ,   0
        delta(i,12)+delta(i,8)*delta(i,11) , delta(i,13)];
    
    
    f= (A*P+P*A'+B*W'+W*B'+2*alpha*P); %  quadratic performance criterion
    %G=10^3-(A*P+P*A'+B*Y'+Y*B'+2*alpha*P);
    g(i,1)=max(eig(f));
    %g(i,2)=max(eig(P));
end
end


 
