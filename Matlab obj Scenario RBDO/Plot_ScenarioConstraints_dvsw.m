function   Plot_ScenarioConstraints_dvsw(obj,d,deltaIdx)
%% DESCRIPTION HERE
%%

Nsep=10;
TheDesing=d(1:obj.Nd);
dLB=TheDesing-abs(TheDesing)*0.1;% plot in the neighbour of d
dUB=TheDesing+abs(TheDesing)*0.1;% plot in the neighbour of d
for i=1:obj.Nd
    d= TheDesing; % back to the original desing
    Di=linspace(dLB(i),dUB(i),Nsep);
    for di=1:length(Di)
        d(i)=Di(di);  % explore design variable di
        Rel=obj.Compute_ReliabilityMetrics(d);
        W{i}(:,di)=Rel.W;
    end
end

for i=1:obj.Nd
    % title('W(\delta,d)')
    subplot(ceil(obj.Nd/3),3,i)
    plot(linspace(dLB(i),dUB(i),Nsep),W{1, i}(deltaIdx,:)','LineWidth',1);
    hold on; grid on;
    xlabel(['d_' num2str(i)]); ylabel('max(g)=w(\delta^(^i^),d)')
    h =line([TheDesing(i) TheDesing(i)], [min(min(W{1, i})) max(max(W{1, i}))]);
    h.LineStyle = '--';
    h.LineWidth=2;
    set(gca,'FontSize',14);
end
end