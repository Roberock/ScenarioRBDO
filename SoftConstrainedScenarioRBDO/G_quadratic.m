function G=G_quadratic(d,delta)
N=size(delta,1); 
nd=size(d,2);  

for i=1:N 
   A_i=hankel(delta(i,1:nd)); 
   B_i=hankel(delta(i,nd+1:2*nd)) ;
   C_i=delta(i,:);
   TempG=d(1)*A_i*d'.^2-d(2)*B_i'*d'+ sin(C_i);
   G(i,:)=TempG(:); 
end

end