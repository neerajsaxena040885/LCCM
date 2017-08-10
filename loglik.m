% Calculates log-likelihood function value for the Latent Class Model
%
% This code is input to Matlab's fminunc command
%
% Input param is a row vector of parameters, dimension 1X[(NF+NV)*NOC+NOC]
%     containing the fixed coefficients, the random coefficient sigma for
%     each latent class and the constant term for each latent class
% Output ll is the scalar value of the negative of the simulated log-likelihood 
%     at the input parameters
% 
%   PURPOSE of this code is to split the param array into fixed beta,
%   random beta and random Ws, which then go into another function llgrad2
%

function [ll, g] =loglik(param)     % The contents of param array are --> {Fixed, Non-standard normal betas (Nv minus std normal) and random Ws (Nv)}. For eg... In case of 1 Fixed and 3 Random, param = {F,B,B,B,W,W,W}
                                    
global NV NF NLCV NOC NAT 

% Splitting the parameter array into fixed, random and latent class
% parameters
for class = 1:NOC
   if NF>0
       % Set of parameter values belonging to a given latent class
       fixvector(class,:) = param(1,(class-1)*(NAT+NV)+1:class*NAT+(class-1)*NV);   % size of fixvector is NOC X NAT
   else
       fixvector(class,:)=[];
   end
   
   if NV>0
       % Set of parameter values belonging to a given latent class
       randvector(class,:) = param(1,class*(NAT+1):class*(NAT+NV));             % size of randvector is NOC X NV
   else
       randvector(class,:)=[];
   end
   
   % Latent class prevalence parameters
   if class < NOC   % Keeping the last latent class as refrence
      lcvector(class,:) = param(1,NOC*(NAT+NV)+NLCV*(class-1)+1:NOC*(NAT+NV)+NLCV*class);   % size of lcvector is NOC X NLCV
   end
end

[p g]=llgrad2(fixvector,randvector,lcvector);   % Giving fixed, random and latent class parameters as input arguments to another function llgrad2

ll=-sum(log(p),2);  % Sum of log-probabilities of all individuals... Extra -ve to make it into convex function for fminunc in Matlab
g=-sum(g,2);        % Sum of gradient of all individuals... Extra -ve to make it into convex function for fminunc in Matlab

%
%

