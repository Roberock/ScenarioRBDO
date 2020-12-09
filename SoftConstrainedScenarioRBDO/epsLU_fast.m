function [epsL, epsU] = epsLU_fast(k,N,bet)
%% Reference Article 
% Title Risk and complexity in scenario optimization 
% Authors S. Garatti and ·M. C. Campi
% Journal Mathematical Programming 
% DOI https://doi.org/10.1007/s10107-019-01446-4
 
% This function provide the lower and upper reliability parameter for a
% convex program as defined in Eq (14) of the refernced paper

% N= number of samples
% k = Number of scenarios \delta_i for which \zeta_i \geq 0, i.e.  f(x,\delta_i) \geq 0
% beta = confidence parameter (e.g. very high confidence beta=10^-8)
alphaL = betaincinv(bet,k,N-k+1); 
alphaU = 1-betaincinv(bet,N-k+1,k);
m1 = [k:1:N];
%aux1 = sum(triu(log(ones(N-k+1,1)*m1),1),2); 
Temp=log(ones(N-k+1,1).*m1');Temp=sort(Temp,'descend');  
Temp(end)=0; CUMSUM=cumsum(Temp); CUMSUM(end)=0;
aux1 =sort(CUMSUM,'descend'); 

%aux2 = sum(triu(log(ones(N-k+1,1)*(m1-k)),1),2); 
Temp=log(ones(N-k+1,1).*(m1-k)');Temp=sort(Temp,'descend');  
Temp(end)=0; CUMSUM=cumsum(Temp); CUMSUM(end)=0;
aux2 =sort(CUMSUM,'descend'); 

coeffs1 = aux2-aux1; m2 = [N+1:1:4*N]; 
%aux3 = sum(tril(log(ones(3*N,1)*m2)),2);
Temp=log(ones(3*N,1).*m2');
Temp=sort(Temp);  CUMSUM=cumsum(Temp); 
aux3 =sort(CUMSUM); 

%aux4 = sum(tril(log(ones(3*N,1)*(m2-k))),2);
Temp=log(ones(3*N,1).*(m2-k)');
Temp=sort(Temp); CUMSUM=cumsum(Temp); 
aux4 =sort(CUMSUM); 

coeffs2 = aux3-aux4;
t1 = 1-alphaL;
t2 = 1;
poly1 = 1+bet/(2*N)-bet/(2*N)*sum(exp(coeffs1 - (N-m1')*log(t1)))... 
    -bet/(6*N)*sum(exp(coeffs2 + (m2'-N)*log(t1)));
poly2 = 1+bet/(2*N)-bet/(2*N)*sum(exp(coeffs1 - (N-m1')*log(t2)))...
    -bet/(6*N)*sum(exp(coeffs2 + (m2'-N)*log(t2)));

if ((poly1*poly2) > 0) 
    epsL = 0;
else    
    while t2-t1 > 1e-10
        t = (t1+t2)/2;
        polyt = 1+bet/(2*N)-bet/(2*N)*sum(exp(coeffs1 - (N-m1')*log(t)))...
            -bet/(6*N)*sum(exp(coeffs2 + (m2'-N)*log(t))); 
        if polyt > 0 
            t1=t;
        else
            t2=t;
        end
    end
    epsL = 1-t2; 
end
t1 = 0; t2 = 1-alphaU;
poly1 = 1+bet/(2*N)-bet/(2*N)*sum(exp(coeffs1 - (N-m1')*log(t1)))...
-bet/(6*N)*sum(exp(coeffs2 + (m2'-N)*log(t1))); 

poly2 = 1+bet/(2*N)-bet/(2*N)*sum(exp(coeffs1 - (N-m1')*log(t2)))... 
    -bet/(6*N)*sum(exp(coeffs2 + (m2'-N)*log(t2))); 

if ((poly1*poly2) > 0) 
    epsL = 0; 
else
    while t2-t1 > 1e-10 
        t = (t1+t2)/2; 
        polyt =1+bet/(2*N)-bet/(2*N)*sum(exp(coeffs1-(N-m1')*log(t)))...
            -bet/(6*N)*sum(exp(coeffs2 + (m2'-N)*log(t)));
        if polyt > 0 
            t2=t; 
        else
            t1=t;
        end
    end
    epsU = 1-t1; 
end
end
