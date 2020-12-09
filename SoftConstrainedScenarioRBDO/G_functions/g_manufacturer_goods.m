function g=g_manufacturer_goods(d,delta,Nk,Atreshold)


% delta = uncertainty in the production of goods
% n = different production resources employed in the production
% q(k,j) =  quantity of resource k used in workplace j to produce a unitary amount of goods 
% d = quantity of goods to be produced 
% Nk = number of resources
% Nd = number of different workplaces where the good is produced  
Nd=length(d);   
g= zeros(size(delta,1),Nk);
for k=1:Nk 
 Qval=1./(delta(:,k).^0.25 +delta(:,Nk+1+(k-1)*Nd:Nk+Nd+(k-1)*Nd));
 g(:,k)=sum(Qval.*d, 2)-Atreshold(k); 
end


% we inform the reader about the mechanism by which qj,k(?) were generated. 
% Let ? = (?1,?2,?1,1,...,?50,1,?1,2,...,?50,2), where, for k = 1,2 and j = 1,...,50, 
% ?k ? U[10,50] (i.e., ?k is uniformly distributed in[10,50]),
% ?j,k ?U[?2.5,2.5], and all these variables are independent
% Then, q_{j,k}(?) = (?_k^{1/4} +?_{j,k})^{?1}, j =1,...,50, k =1,2. 
% Moreover, in the simulation, we took a1 = a2 = 1

end