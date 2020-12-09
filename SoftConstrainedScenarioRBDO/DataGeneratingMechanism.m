function delta=DataGeneratingMechanism(N)
% Composed by inliers bi-variate Gaussian with different % of outliers

% N = number of samples e.g. N=1e5;
% delta= vector of uncertain scenarios
% N= number of samples

%% START
Corr=0.5; % correlation values
MeanDelta=[-2,-2]; % mean values

InlierPercent=0.97; % 97% inliers and 3% are outlier
rng default; % reset random seed
%Count=0;
delta_inliers= mvnrnd(MeanDelta,[1 Corr; Corr, 2^2],floor(N*InlierPercent)) ;
Count=size(delta_inliers,1);
delta_outliers=[];
while Count<N % generate samples until N are obtained
    r=rand(); 
        if r>0.80% 20 % of the outliers are outliers of type 1
            delta_outliers=[delta_outliers;mvnrnd([-5 -5],[2 0.4; 0.4, 2],[1])]; % Generate outliers type 1
        else  %  the remaining percentage of outliers are outliers of type 2
            delt1=unifrnd(normrnd(-2,0.5),normrnd(3,0.5));
            delta_outliers=[delta_outliers;[delt1 +(abs(delt1)*rand())^2*(1+rand())]];  % Generate outliers type 2
        end  
    Count=size(delta_inliers,1)+size(delta_outliers,1);
end
delta=[delta_outliers;delta_inliers];
delta=delta(1:N,:);
end