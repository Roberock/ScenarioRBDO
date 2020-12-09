function Area=ApproxFeasibilityArea(delta,alpha,lambda,g_fun,Delta_idx,LBd,UBd,Ngridpoints)
%% Numerical approximation of the area s.t. Var_alpha<lambda and CVar_alpha<lambda for individual requirements
%Ngridpoints=25;
D1vec=linspace(LBd(1) , UBd(1) ,Ngridpoints);
D2vec=linspace(LBd(2) ,UBd(2) ,Ngridpoints);
[D1,D2]=meshgrid(D1vec,D2vec);
Designs= [D1(:) D2(:)]; % vector of D coordinates

%% Prealocate
if ~isempty(Delta_idx)
    [VAR_alpha_g1,VAR_alpha_g2,CVAR_alpha_g1,CVAR_alpha_g2]=deal(zeros(size(Designs,1),length(Delta_idx) ));
    [VAR_alpha_w,CVAR_alpha_w ]=deal(zeros(size(Designs,1),length(Delta_idx)));
else
    [VAR_alpha_g1,VAR_alpha_g2,CVAR_alpha_g1,CVAR_alpha_g2]=deal(zeros(size(Designs,1),1));
    [VAR_alpha_w,CVAR_alpha_w ]=deal(zeros(size(Designs,1),1));
end
%% Analyse the desings and the removed samples
for i=1:size(Designs,1)
    G=g_fun(Designs(i,:),delta);
    W=max(G,[],2);
    if ~isempty(Delta_idx)
        for k=1:length(Delta_idx)
            % joint req
            VAR_alpha_w(i,k) =prctile(W([1:Delta_idx(k)-1,Delta_idx(k)+1:end]),alpha*100);
            CVAR_alpha_w(i,k) =CVaR(W([1:Delta_idx(k)-1,Delta_idx(k)+1:end]),alpha);
            % individual req
            VAR_alpha_g1(i,k)=prctile(G([1:Delta_idx(k)-1,Delta_idx(k)+1:end],1),alpha *100);
            VAR_alpha_g2(i,k)=prctile(G([1:Delta_idx(k)-1,Delta_idx(k)+1:end],2),alpha *100);
            CVAR_alpha_g1(i,k)=CVaR(G([1:Delta_idx(k)-1,Delta_idx(k)+1:end],1),alpha );
            CVAR_alpha_g2(i,k)=CVaR(G([1:Delta_idx(k)-1,Delta_idx(k)+1:end],2),alpha );
        end
    else
        % joint req
        VAR_alpha_w(i) =prctile(W,alpha*100);
        CVAR_alpha_w(i) =CVaR(W,alpha);
        % individual req
        VAR_alpha_g1(i)=prctile(G(:,1),alpha*100); 
        VAR_alpha_g2(i)=prctile(G(:,2),alpha*100);
        CVAR_alpha_g1(i)=CVaR(G(:,1),alpha);
        CVAR_alpha_g2(i)=CVaR(G(:,2),alpha);
        
        IsFeasible_var_g1(i)= VAR_alpha_g1(i)<0;
        IsFeasible_var_g2(i)= VAR_alpha_g2(i)<0;
        IsFeasible_cvar_g1(i)= CVAR_alpha_g1(i)<0;
        IsFeasible_cvar_g2(i)= CVAR_alpha_g2(i)<0;
    end
end
% collect results
Area.VarW=mean(VAR_alpha_w<lambda(1) );
Area.CVarW=mean(CVAR_alpha_w<lambda(1) );
Area.Varg(1,:)=mean(VAR_alpha_g1<lambda(1));
Area.Varg(2,:)=mean(VAR_alpha_g2<lambda(2));
Area.CVarg(1,:)=mean(CVAR_alpha_g1<lambda(1));
Area.CVarg(2,:)=mean(CVAR_alpha_g2<lambda(2));
Area.Contour.D1vec=D1vec;
Area.Contour.D2vec=D2vec;
Area.Contour.VAR_alpha_w=reshape(VAR_alpha_w,Ngridpoints,Ngridpoints);
Area.Contour.CVAR_alpha_w=reshape(CVAR_alpha_w,Ngridpoints,Ngridpoints);
Area.Contour.VAR_alpha_g1=reshape(VAR_alpha_g1,Ngridpoints,Ngridpoints);
Area.Contour.VAR_alpha_g2=reshape(VAR_alpha_g2,Ngridpoints,Ngridpoints);
Area.Contour.CVAR_alpha_g1=reshape(CVAR_alpha_g1,Ngridpoints,Ngridpoints);
Area.Contour.CVAR_alpha_g2=reshape(CVAR_alpha_g2,Ngridpoints,Ngridpoints); 

Area.Contour.IsFeasible_var_g1=reshape(IsFeasible_var_g1,Ngridpoints,Ngridpoints);
Area.Contour.IsFeasible_var_g2=reshape(IsFeasible_var_g2,Ngridpoints,Ngridpoints); 
 
Area.Contour.IsFeasible_cvar_g1=reshape(IsFeasible_cvar_g1,Ngridpoints,Ngridpoints);
Area.Contour.IsFeasible_cvar_g2=reshape(IsFeasible_cvar_g2,Ngridpoints,Ngridpoints); 
 
end