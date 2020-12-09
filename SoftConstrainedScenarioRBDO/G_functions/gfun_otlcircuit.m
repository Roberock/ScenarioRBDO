function [g] = gfun_otlcircuit(d,delta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OTL CIRCUIT FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RANGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rb1 ? [50, 150]	resistance b1 (K-Ohms)
% Rb2 ? [25, 70]	resistance b2 (K-Ohms)
% Rf ? [0.5, 3]     resistance f  (K-Ohms)
% Rc1 ? [1.2, 2.5]	resistance c1 (K-Ohms)
% Rc2 ? [0.25, 1.2] resistance c2 (K-Ohms)
% beta ? [50, 300]	current gain  (Amperes)
%%%OUTPUT AND INPUT:%%%%%%%%%%%%%%%%%%%%%%%%vvv
 
% Vm = midpoint voltage
% d = [Rb1, Rb2, Rf, Rc1, Rc2, beta]
% delta = [errRb1, errRb2, errRf, errRc1, errRc2, errbeta] errors on the desing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rb1  = d(1)+delta(:,1);
Rb2  = d(2)+delta(:,2);
Rf   = d(3)+delta(:,3);
Rc1  = d(4)+delta(:,4);
Rc2  = d(5)+delta(:,5);
beta = d(6)+delta(:,6);

Vb1 = 12*Rb2 ./ (Rb1+Rb2);
term1a = (Vb1+0.74) .* beta .* (Rc2+9);
term1b = beta.*(Rc2+9) + Rf;
term1 = term1a ./ term1b;

term2a = 11.35 .* Rf;
term2b = beta.*(Rc2+9) + Rf;
term2 = term2a ./ term2b;

term3a = 0.74 .* Rf .* beta .* (Rc2+9);
term3b = (beta.*(Rc2+9)+Rf) .* Rc1;
term3 = term3a ./ term3b;

Vm = term1 + term2 + term3;

g1 = term1-4.5;
g2 =term2-0.01;
g3 = term3-0.5;
g4= Vm-6;
g = [g1,g2,g3,g4];
end
