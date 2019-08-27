%% Basic GA parameters
gaDat.Objfun='objfun_schwefel';
lb=[-500 -500];
ub=[500 500];
gaDat.FieldD=[lb; ub];
% Execute GA
gaDat=ga(gaDat);
% Result are in
gaDat.xmin
gaDat.fxmin
