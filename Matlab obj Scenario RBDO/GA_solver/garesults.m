function garesults(gaDat)
% Optional user task executed when the algorithm ends

% For instance, final result
disp('------------------------------------------------')
disp('######   RESULT   #########')
disp(['   Objective function for xmin: ' num2str(gaDat.fxmin)])
disp(['   xmin: ' mat2str(gaDat.xmin)])
disp('------------------------------------------------')