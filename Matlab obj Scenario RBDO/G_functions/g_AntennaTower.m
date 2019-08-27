function g=g_AntennaTower(delta,d)
%% Antenna Tower (Broggi Patelli etc.) see OpenCossan Documentation 
% a robust optimization of a 25-bars truss structure,% i.e. an antenna tower, is carried out. The direction of the force acting
% on the  structure and the structural parameters are affected by% uncertainties.
% The objective of the optimization is to minimize the volume of the
% structure (proportional to the costs of the structure), while assuring
% that the maximum nodal displacement is under a certain threshold. 
 
%% There are 28 iid random quantities delta representing :
% Young's moduli of the beams i=1,...,25

% Create a set of 2 uniform distributed random variables for the force
% direction. This direction is a spherical angle deviation of +- 5 degrees
% from the vertical direction, and a totally random direction in the
% horizontal plane.
 
Emodulus=delta(:,1:25);
phi=delta(:,26);
theta=delta(:,27);
%% This design variables are 6 groups of identical beams characterized by the same section:
% d1 = Parameter('value',0.4);
% d2 = Parameter('value',0.1);
% d3 = Parameter('value',3.4);
% d4 = Parameter('value',1.3);
% d5 = Parameter('value',0.9);
% d6 = Parameter('value',1.0)
 
%% TRUSSMAXDISPSCRITP - Main MIO code, script version
% for each sample
FX=-100e3.*cos(phi).*sin(theta);
FY=-100e3*sin(phi).*sin(theta);
FZ=-100e3.*cos(theta); 

for sam = 1:size(delta,1)
    %% construct the antenna structure from the sampled values. D is a
    % structure containing:
    % - nodal coordinate (D.Coord)
    % - nodal connectivity (D.Con)
    % - boundary conditions (D.Re)
    % - nodal loads (D.Load)
    % - beams Young's moduli (D.E)
    % - beams sections (D.A)
    %%  Definition of the 25 beams truss structure
    % Get structural parameter from input struct sample
    % E - random variables (a for-loop with eval is use here for brevity)
 
    % Fx, Fy, Fz - function of the random variables theta and phi
    Fx = FX(sam,:);    Fy = FY(sam,:);    Fz = FZ(sam,:);
    % Ai (sections of the beams) - parameter. The values of these parameters
    % are set from the outer loop design variables
    A1 = d(1);    A2 = d(2);
    A3 = d(3);    A4 = d(4);
    A5 = d(5);    A6 = d(6);
    
    %% Construct the struct of the truss
    %  Nodal Coordinates
    Coord=[-37.5 0 200;37.5 0 200;-37.5 37.5 100;37.5 37.5 100;37.5 -37.5 100;...
        -37.5 -37.5 100;-100 100 0;100 100 0;100 -100 0;-100 -100 0];
    
    %  Connectivity
    Con=[1 2;1 4;2 3;1 5;2 6;2 4;2 5;1 3;1 6;3 6;4 5;3 4;5 6;...
        3 10;6 7;4 9;5 8;4 7;3 8;5 10;6 9;6 10;3 7;4 8;5 9];
    
    % Definition of Degree of freedom (free=0 &  fixed=1).
    % The nodes 7 to 10 are fixed in all the three degrees of freedom
    Re=zeros(size(Coord));Re(7:10,:)=[1 1 1;1 1 1;1 1 1;1 1 1];
    
    % Definition of Nodal loads. The Loads are applied at the two top nodes.
    Load=zeros(size(Coord));Load(1:2,:)=[Fx Fy Fz; Fx Fy Fz];
    
    % Definition of Modulus of Elasticity
    E = Emodulus(sam,:);
    
    % Definition of beams sections
    A=[A1 A2 A2 A2 A2 A3 A3 A3 A3 A1 A1 A4 A4 A5 A5 A5 A5 A6 A6 A6 A6 A3 A3 A3 A3];
    
    % Convert to structure
    D=struct('Coord',Coord','Con',Con','Re',Re','Load',Load','E',E','A',A');
    %% compute the nodal displacements in the 3 degrees of freedom
    %   Compute the stresses, nodal diplacements and reaction forces given the
    %   struct describing the truss. 
    %   The stiffness matrix is assembled and used to computed the
    %   aforementioned quantities.
    w=size(D.Re);S=zeros(3*w(2));U=1-D.Re;f=find(U);
    for i=1:size(D.Con,2)
        H=D.Con(:,i);C=D.Coord(:,H(2))-D.Coord(:,H(1));Le=norm(C);
        T=C/Le;s=T*T';G=D.E(i)*D.A(i)/Le;Tj(:,i)=G*T;
        e=[3*H(1)-2:3*H(1),3*H(2)-2:3*H(2)];S(e,e)=S(e,e)+G*[s -s;-s s];
    end
    U(f)=S(f,f)\D.Load(f);F=sum(Tj.*(U(:,D.Con(2,:))-U(:,D.Con(1,:))));
    R=reshape(S*U(:),w);R(f)=0;
    %% compute the norms of the nodal displacements
    normU = zeros(1,size(U,2));
    for inode = 1:size(U,2)
        normU(inode) = norm(U(:,inode));
    end
    % assign the maximum displacement to the first output
    maxDisp = max(normU);
    %% G1-6
    g(sam,:)=[normU(1:3)-5.5 normU(4:6)-0.5]; % max displacement =0.3
     
    % compute the beam length
%       Compute length of beams given the struct describing the truss. 
%     beamLengths=zeros(1,size(D.Con,2));
%     for i=1:size(D.Con,2)
%        H=D.Con(:,i);C=D.Coord(:,H(2))-D.Coord(:,H(1));
%        beamLengths(i)=norm(C);
%     end
%     % compute the beam volumes and assign them to the sencond output
%     beamVolumes(sam) = beamLengths' .* D.A;
    
end
 
end
