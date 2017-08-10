% Calculate logit probability for chosen alternatives for each person
%    with multiple choice situations for each person and multiple people,
%    using globals for all inputs except coefficients.
%
% Simulated Mixed Logit Probability and gradient alongwith Latent class
% component
%
% IDEA for the latent class model is to evaluate an error component logit
% model for each set of parameters. An explanation of Logit probability
% given by Prof. Train is given below:
%
% Logit probability is Prod_t { exp(V_*t)/SUM_j[exp(V_jt)] } 
%             = Prod_t { 1 / (1+ Sum_j~=* [exp(V_jt-V-*t)] }
%    where * denotes chosen alternative, j is alternative, t is choice situation.
% Using differences from chosen alternative reduces computation time (with
%    one less alternative), eliminates need to retain and use the dependent
%    variable, and avoids numerical problems when exp(V) is very large or
%    small, since prob is 1/(1+k) which can be evaluated for any k, including
%    infinite k and infinitesimally small k. In contrast. e(V)/sum e(V) can
%    result in numerical "divide by zero" if denominator is sufficiently small
%    and NaN if numerator and denominator are both numerically zero.
% 
% Input fixvector contains a vector of fixed coefficients, and has dimension NOCXNF.
% Input randvector contains a vector of fixed coefficients, and has dimension NOCXNV.
% Input lcvector contains a vector of fixed coefficients, and has dimension NOCX1. 
% Output p contains the logit probabilities, which is a row vector of dimension 1xNP.
% Output g contains the gradients of log(p), which is a matrix
% [NOC*(NF+NV)+NOC] x NP.
%
% These matrices will be then passed back to loglik function which will sum
% up both the matrices along the NP dimension.
% Code assumes that all GLOBALS already exist.


function [p, g]=llgrad2(fixvector,randvector,lcvector); 

global NV NF NP NDRAWS NALTMAX NCSMAX NMEM NOC NLCV 
global X S XF XLCV

p=zeros(1,NP);
g=zeros(NOC*(NF+NV)+(NOC-1)*NLCV,NP);
sumppp=zeros(1,NP,NMEM);

% Calculating the product of beta coefficients (levector) and
% covariates
vlcv=lcvector*XLCV; % matrix multiplication.. lcvector [(NOC-1) X NLCV] and XLCV [NLCV X NP].. Resulting matrix is [(NOC-1) X NP]
vlcv=exp(vlcv);     % taking exponentials
plcv=1./(1+sum(vlcv,1)); % Evaluating the probability of base latent class by summing up the exponentials along the first dimension.. Resulting matrix is 1 X NP

for class=1:NOC   %Finding the likelihood and gradient functions using different latent class parameters 

    % Evaluating the probabilities of the latent class prevalence
    if class < NOC
        PPLCV(class,:)=vlcv(class,:).*plcv;
    else
        PPLCV(class,:)=plcv;   % The last latent class is set to base
    end
    
    v=zeros(NMEM,NALTMAX-1,NCSMAX,NP); % Observed part of utility... Sum of fixed and error component parts
    if NF > 0
        ff=reshape(fixvector(class,:),1,1,NF,1);
        ff=repmat(ff,[NALTMAX-1,NCSMAX,1,NP]);   % Fixed array therefore no NDRAWS
        vf=reshape(sum(XF.*ff,3),NALTMAX-1,NCSMAX,NP);  %vf is (NALTMAX-1) x NCSMAX x NP     % Vf is the observed utility due to fixed coefficients. Vf = Beta*Xf  % XF.*FF is element-element multiplication. 3 is the third dimension of XF for which multiplication happens
    else
        vf=zeros(NALTMAX-1,NCSMAX,NP);
    end
    vf=repmat(vf,[1,1,1,NMEM]);
    
    if NV > 0
        v=(repmat(X,[1,1,1,1,1])*randvector(class,1)); %v is (NALTMAX-1) x NCSMAX x NV x NP x NMEM..... Multiplying the matrix with a scalar
        v=reshape(sum(v,3),NALTMAX-1,NCSMAX,NP,NMEM);             %v is (NALTMAX-1) x NCSMAX x NP x NMEM     % Adding the random terms... This reduces the dimensions to 4.. Doesn't make any difference as there is just one random term.. Now it can be added to Vf
        v=v+vf;      % Total observed utility (fixed + error component) for an alternative, choice situation, individual and draw
    else
        v=vf;
    end
    
    v=exp(v);   % Taking exp of the observed utility
    v(isinf(v))=10.^20;  %As precaution when exp(v) is too large for machine       % Defining Infinity as 10.^20.. Setting an upper bound
    v=v.*repmat(S,[1,1,1,NMEM]);    % This command multiplies total utility with alternative's availability... This command is not important in cases where all alternatives are available everytime.
    eval(['V' num2str(class) '=v;']);   % Capture the value of v for each latent class
    pp=1./(1+sum(v,1)); %pp is 1 x NCSMAX x NP x NMEM   % Gives the probability of making t choices for every person across draws.. For example, for each person and draw, a set of 10 choice probabilities are generated (for the chosen alternative)
    % Capturing pp value for each latent class
    eval(['PP' num2str(class) '=pp;']);
    eval(['PPP' num2str(class) '=reshape(pp,NCSMAX,NP,NMEM);']);    %PPP is NCSMAX x NP x NMEM
    eval(['PPP' num2str(class) '=prod(PPP' num2str(class) ',1);']);    %PPP is now 1xNPxNMEM    % Taking the product of all probabilites across t choice tasks for an individual
    eval(['sumppp=sumppp + PPP' num2str(class) '.*repmat(PPLCV(class,:),[1,1,NMEM]);']);     % Accumulating the product class(i)*p(i)
end
%Back to prob
pp=sum(sumppp,3);       %pp is 1xNP.... Its a summation of class(i)*p(i)
p=p+pp;
p=p./NDRAWS;    % taking average probability over NDRAWS
p(1,isnan(p))=1; %Change missing values to 1, as a precaution.

%-------------------------
%Calculate gradient
%-------------------------

for class=1:NOC
  if class < NOC    % I now have NOC-1 latent class parameters
      I=ones(NLCV,NP,NMEM);
      eval(['gglcv=repmat(PPP' num2str(class) ',[NLCV,1,1]).*repmat(XLCV,[1,1,NMEM]);']);
      prodpplcv=repmat(PPLCV(class,:),[NLCV,1,NMEM]).*(I-repmat(PPLCV(class,:),[NLCV,1,NMEM]));
      eval(['grlc' num2str(class) '=(gglcv.*prodpplcv)./repmat(sumppp,[NLCV,1,1]);']);  % Size of the matrix is 1 x NP x NMEM. Calculating gradient wrt to latent class. Storing each of them separately. 
  end
  
  % Preparing for the fixed part gradient calculation now...
  eval(['gg=V' num2str(class) '.*repmat(PP' num2str(class) ',[NALTMAX-1,1,1,1]);']);    %Probs for all nonchosen alts NALTMAX-1 x NCSMAX x NP x NMEM    % AWESOME: This is not the gradient... It's just a probability
  % gg=v.*repmat(pp,[NALTMAX-1,1,1,1]);   
  gg=reshape(gg,NALTMAX-1,NCSMAX,1,NP,NMEM);
  if NF>0
      grf=-repmat(gg,[1,1,NF,1,1]).*repmat(XF,[1,1,1,1,NMEM]);  % this is correct.. See the proof in my notebook... This derivative represents wrt fixed coefficient
      grf=reshape(sum(sum(grf,1),2),NF,NP,NMEM);   %NFxNPxNMEM  % Total gradient for a fixed parameter is [sum over alternatives {sum over choice tasks (XF * Tot prob v)}].. Sum of 2*10 probability values is one cell of this array.
      eval(['grfix' num2str(class) '=(repmat(PPP' num2str(class) ',[NF,1,1]).*grf.*repmat(PPLCV(class,:),[NF,1,NMEM]))./repmat(sumppp,[NF,1,1]);']);
  else
      eval(['grfix' num2str(class) '=[];']);
  end
   
  if NV>0
      grw=-repmat(gg,[1,1,NV,1,1]).*X;
      grw=reshape(sum(sum(grw,1),2),NV,NP,NMEM);    %NFxNPxNMEM
      eval(['grrand' num2str(class) '=(repmat(PPP' num2str(class) ',[NV,1,1]).*grw.*repmat(PPLCV(class,:),[NV,1,NMEM]))./repmat(sumppp,[NV,1,1]);']);
  else
      eval(['grfix' num2str(class) '=[];']);
  end
end

% Preparing the final matrix of gradients of size NOC*(NF+NV)+(NOC-1)*NLCV X NP X
% NMEM
for class=1:NOC
    if class == 1
        eval(['gr=[grfix' num2str(class) ';grrand' num2str(class) '];']);
        continue;
    end
    eval(['gr=[gr;grfix' num2str(class) ';grrand' num2str(class) '];']);    % Continue appending gradients into gr matrix
end
for class=1:(NOC-1)     % Only NOC-1 latent class parameters... The last latent class is the base
   eval(['gr=[gr;grlc' num2str(class) '];']); % Appending the latent class prevalence gradient
end


%Gradient
   % CHECK WITH DUBEYJI... I think the code below need not be activated...
   % See my derivations!!
   %gr=gr.*repmat(pp,[NOC*(NF+NV)+NOC,1,1]);    % Populating a single big matrix of gradient for all parameters put together. The number of paramters are --> NOC*(NF+NV)+NOC
   g=g+sum(gr,3);  %gr is (NF+NV+NV) x NP.... The last  dimension NMEM went away by adding it up
 
%Gradient
   g=g./NDRAWS; % taking average probability over NDRAWS
   %g=g./repmat(p,NOC*(NF+NV)+NOC,1);   % Dividing gradient by probability again.. Don't know.. Why this is happening?? ASK DUBEYJI AGAIN!!!!!!