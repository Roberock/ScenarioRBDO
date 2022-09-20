
%% Example of how to use Conditional-Value at Risk (CVaR) 
MAtrixVoltageDeviations= abs(normrnd(0,0.01,[14,500])'); % your matrix of samples e.g. voltage deviations from reference
AbsVoltDev=abs(MAtrixVoltageDeviations); % take the absolute value
VectorAbsVoltDev=AbsVoltDev(:); % the rows vector of voltage deviations [[Number of busses * Number of samples] x 1]

%% Now compute the CVaR
Samples=VectorAbsVoltDev; % input samples 
Alpha=0.8;  % select an alpha level alpha in [0,1]  
[CVAR_alpha_80]=CVaR(Samples,Alpha);
% NOTE-1: The valeu 100*(1-alpha) gives you the percentage of hihg-voltage sceanrio profile you seek to optimize
% NOTE-2: A value alpha=0 should be equal to the mean Value of the
% distribution (apart from minor numerical approx error). 
disp(['is the value CVaR(alpha=0)=' num2str(CVaR(Samples,0)) ', equal to the mean voltage deviation mu=' num2str(mean(VectorAbsVoltDev)), ' ?'])
%% Let's have a closer look to the Value at Risk (The CVaR function starts by doing this)
[F,X]=ecdf(Samples); % get emirical comulative distribution function  
i=find(F >= Alpha  , 1, 'first'); % find the first value of the samples corresponding to the alpha-level
XVaR=X(i); % the Value at Risk (the Voltage deviation corresponding to the alpha-percentile)
%% Visulaize the difference between VaR and CVaR Plot
subplot(2,1,1)

ecdf(VectorAbsVoltDev); xlabel('Absolute Voltage deviation from 1 p.u.'); grid on;
xline(XVaR,'-.b') ;title('Comulative Distribution Functin')
xline(CVAR_alpha_80,'--r') 
legend('Empirical CDF','Value at Risk', 'CVaR')
subplot(2,1,2)

ksdensity(VectorAbsVoltDev);xlabel('Absolute Voltage deviation from 1 p.u.'); grid on;
xline(XVaR,'-.b') ;title('Kernel density estimator')
xline(CVAR_alpha_80,'--r') 
legend('PDF','Value at Risk', 'CVaR')



