function gaiteracion(gaDat)
%  Optional user task executed at the end of each iteration
%

% For instance, results of the iteration
 disp('------------------------------------------------')
 disp(['Iteration: ' num2str(gaDat.gen)])
 disp(['   xmin: ' mat2str(gaDat.xmin) ' -- f(xmin): ',num2str(gaDat.fxmin)])
%% --------------------------------------------------------
