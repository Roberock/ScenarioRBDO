function [A,B,Aeq,Beq, LB, UB,options]=getFminconSet(obj)
% load the optimization structures from the scenario RBDO object
UB=obj.OptimizerData.UB;
LB=obj.OptimizerData.LB;
Beq=obj.OptimizerData.Beq;
Aeq=obj.OptimizerData.Aeq;
B=obj.OptimizerData.B;
A=obj.OptimizerData.A; 
options=obj.OptimizerData.options;  
end
