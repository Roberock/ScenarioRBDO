function [g_p]=g_ControllerCrespo(delta,d)
%% Input
% delta      = the available scenarios/samples [N samples x N uncertain factors]
% d     = the design of the controller [1x N desvar]
% do_linear  = the model to be employed 1=low fidelity linear 0=high fidelity non linear
% tvec       = time steps vector
%% Output
% g_p    = performance functions evaluated in (delta,desing) 
% [N samples x N performance functions]  g>0 is failure g<=0 is safe  
%%

do_linear=1; 
tvec=0:0.002:16; 
show=0; % do not plot
%[I_f]=zeros(1,size(psams,1));
 
    parfor i=1:size(delta,1)
        g_p(i,:)=evaluate_g2(delta(i,:),d,do_linear,tvec,show); 
    end 
  
end  