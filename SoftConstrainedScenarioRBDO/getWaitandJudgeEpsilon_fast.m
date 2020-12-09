function out = getWaitandJudgeEpsilon_fast(SN,N,bet)

% SN = number of support scenarios
% N = number of scenarios available
% bet= confidence level

% out= epsilon(1), epsilpon(2),.....,epsilon(SN) % vector of bounds for
% different number of supports k=1,...,SN

out = zeros(SN+1,1);
for k = 0:SN
    m = [k:1:N];
    Temp = sort(log(ones(N-k,1).*m(2:end)'),'descend');
    aux1 = [sort(cumsum(Temp),'descend');0];
    Temp = sort(log(ones(N-k,1).*(m(2:end)-k)'),'descend');
    aux2 = [sort(cumsum(Temp),'descend');0]; 
    coeffs = aux2-aux1;
    t1 = 0;
    t2 = 1;
    while t2-t1 > 1e-10
        t = (t1+t2)/2;
        val = 1 - bet/(N+1)*sum( exp(coeffs-(N-m')*log(t)) );
        if val >= 0
            t2 = t;
        else
            t1 = t;
        end
    end
    out(k+1) = 1-t1;
end

 
