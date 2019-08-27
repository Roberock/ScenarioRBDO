function plot_2D_SafeFailDomain(obj,d,DeltaIndex,GIndex,PlotGall)
%% DESCRIPTION HERE 
%%
p_reduced=obj.delta;
plims=[min(p_reduced)-abs(min(p_reduced)).*0.1;max(p_reduced)+abs(max(p_reduced)).*0.1]';  % define upper and lower truncation on the delta space
%% Show 2d-sections of individual and overall safe/failure domains and superimpose observations
% activepos=[1,2]; %components of p being varied. Choose any 2 integers between 1 and 6
nsteps=100;
vec1=plims(DeltaIndex(1),1):(plims(DeltaIndex(1),2)-plims(DeltaIndex(1),1))/nsteps:plims(DeltaIndex(1),2);
vec2=plims(DeltaIndex(2),1):(plims(DeltaIndex(2),2)-plims(DeltaIndex(2),1))/nsteps:plims(DeltaIndex(2),2); 
g_i=cell(1,obj.Ng);
for i=1:length(vec1)
    for j=1:length(vec2)
        p=p_reduced(1,:);
        p(DeltaIndex(1))=vec1(i); % assign the i-th value to factor activepos(1)
        p(DeltaIndex(2))=vec2(j); % assign the j-th value to factor activepos(2)
        g=obj.g_fun(p,d(1:obj.Nd));%g=evaluate_g function in p and for d
        for ng=1:obj.Ng
            g_i{ng}(j,i)=g(ng);
        end
        g_all(j,i)=max(g);
    end
end
if  PlotGall==1
  % plot overall G
ax=figure(1);
hold on
if max(max(g_all))<0 % no failure
    surf(vec1,vec2,sign(g_all));shading flat; colormap(ax,[ 0 0  0.5625; 0.5  0 0 ]); view(2);
elseif min(min(g_all))>0 % all failurs
    surf(vec1,vec2,sign(g_all));shading flat;colormap(ax,[ 0 0  0.5625; 0.5  0 0 ]) ;view(2);
else
    surf(vec1,vec2,sign(g_all));shading flat;colormap(ax,[ 0 0  0.5625; 0.5  0 0 ])  ;view(2);
end
w=scatter3(p_reduced(:,DeltaIndex(1)),p_reduced(:,DeltaIndex(2)),p_reduced(:,DeltaIndex(2))*0+2,'.w'); 
w.MarkerEdgeAlpha =0.5;
axis tight; view(2);
xlabel(['\delta_' num2str(DeltaIndex(1))]);
ylabel(['\delta_' num2str(DeltaIndex(2))]);
title('All Reqs');
set(gca,'FontName','Arial','FontSize',20);  
    
elseif PlotGall==0
figure
for ng=GIndex % plot single Failure regions
    ax=subplot(2,ceil(length(GIndex)/2)+1,ng);
    if max(max(g_i{ng}))<0 % no failure
        surf(vec1,vec2,sign(g_i{ng}));shading flat; colormap(ax,[ 0 0  0.5625; 0.5  0 0 ]); view(2);
    elseif min(min(g_i{ng}))>0 % all failurs
        surf(vec1,vec2,sign(g_i{ng}));shading flat;colormap(ax,[ 0 0  0.5625; 0.5  0 0 ]);view(2);
    else
        surf(vec1,vec2,sign(g_i{ng}));shading flat;colormap(ax,[ 0 0  0.5625; 0.5  0 0 ]);view(2);
    end
    xlabel(['\delta_' num2str(DeltaIndex(1))]);
    ylabel(['\delta_' num2str(DeltaIndex(2))]);title(['g_' num2str(ng)]);
    set(gca,'FontName','Arial','FontSize',20);
    axis tight;
end
% plot overall G
ax=subplot(2,ceil(length(GIndex)/2)+1,2*ceil(length(GIndex)/2)+2);
hold on
if max(max(g_all))<0 % no failure
    surf(vec1,vec2,sign(g_all));shading flat;colormap(ax,[ 0 0  0.5625; 0.5  0 0 ]); view(2);
elseif min(min(g_all))>0 % all failurs
    surf(vec1,vec2,sign(g_all));shading flat;colormap(ax,[ 0 0  0.5625; 0.5  0 0 ]);view(2);
else
    surf(vec1,vec2,sign(g_all));shading flat;colormap(ax,[ 0 0  0.5625; 0.5  0 0 ]);
end
w=scatter3(p_reduced(:,DeltaIndex(1)),p_reduced(:,DeltaIndex(2)),p_reduced(:,DeltaIndex(2))*0+2,'.w'); 
w.MarkerEdgeAlpha =0.5;
axis tight; view(2);
xlabel(['\delta_' num2str(DeltaIndex(1))]);
ylabel(['\delta_' num2str(DeltaIndex(2))]);
title('All Reqs');
set(gca,'FontName','Arial','FontSize',20);
end

end
