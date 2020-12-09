function g=G_fun_Shueller2004(d,delta)

%% see also Rocchetta at al. https://doi.org/10.1016/j.ress.2019.106755
% input
% d = desing vector
% delta = matrix of uncertain parameters  [Number of samples x Number of uncertain factors]
% output
% g = realizations of the reliability performance function  
  
%% G1: 
 g(:,1)= d(1)+d(2).^4.*(delta(:,1)-delta(:,2)).^4-(delta(:,1)-delta(:,2))./sqrt(2) ;
  
%% G2: 
 g(:,2)= d(1)+d(2).^4.*(delta(:,1)-delta(:,2)).^4+(delta(:,1)-delta(:,2))./sqrt(2);
  
%% G3: 
 g(:,3)= d(3).*(d(4).*delta(:,1)-delta(:,2))+(d(2)*5.682)./sqrt(2)+2.2;
  
%% G4: 
 g(:,4)= d(3).*(delta(:,2)-d(4).*delta(:,1))+(d(2)*5.682)./sqrt(2)+2.2;
 
 g=-g; % a failure if g>0
end