function Plot_Bounds_and_groundtruth(Epsilonbounds,X,Violation_groundtruth)

%Epsilonbounds= lower and upper bounds on the true violation probability
%Violation_groundtruth = true probability of violation (e.g estimated via Monte Carlo)

%%
patch([X fliplr(X)], [Epsilonbounds(1,1:length(X)) fliplr(Epsilonbounds(2,1:length(X)))], 'r','Facealpha',0.1)
grid on
box on
xlabel('Number of samples, N')
ylabel('Violation Probability')
hold on
plot(X,Violation_groundtruth(1:length(X)),'k','LineWidth',2)
legend('$[\underline{\epsilon}(k),\overline{\epsilon}(k)]$','$V(d^\star)$','Interpreter','latex')
set(gca,'FontSize',18)
xlim([min(X),max(X)])

end