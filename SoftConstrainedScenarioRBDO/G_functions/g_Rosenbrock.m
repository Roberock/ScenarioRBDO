
function g=g_Rosenbrock(delta,d)
% Rosenbrock FUNCTION
%
% Test Function Series for Optimisation
%
% Copyright (C) 2016 Leong Kuan Yew
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%    * Neither the name of the Volkswagen AG nor the names
%      of its contributors may be used to endorse or promote products derived
%      from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
%
% Domain design: -30 <= xi <= 30
% delta are distributed arround 1 (the optimum in a deterministic Rosenbrock)
% Minimun = 0
% Input:
% delta = [xi...xn], for i = 1, n = 2 by default
%

Ndelta=size(delta,2); % default n=2
Ndesign=length(d); 
g = zeros(size(delta,1),Ndelta-1);
for j=1:Ndesign
    for i=1:Ndelta-1
        g(:,i) = g(:,i)+ (15.*(delta(:,i+1) - d(j).^2).^2 + (delta(:,i) - 1).^2);
    end 
end
g=g-60;% g<50 is the reliability criteria 
end