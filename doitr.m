% Create Global variables to use in estimation
cp=XMAT(:,1); % person
 
nn=zeros(NCS,1);
for n=1:NCS;
    nn(n,1)=sum(XMAT(:,2) == n,1);  %This condition counts the number of lines that have the same choice taskid
end;
NALTMAX=max(nn);    % This gives the number of alternatives 

% Total number of alternatives faced by all people in all choice situations combined.
% This is the number of rows of data in XMAT below.
NROWS=1450.*NALTMAX;

nn=zeros(NP,1);
for n=1:NP;
   k=(XMAT(:,1)==n);
   k=XMAT(k,2);
   nn(n,1)=1+k(end,1)-k(1,1);
end;
NCSMAX=max(nn);  %Maximum number of choice situations faced by any person

% Setting the state of a random number generator. For more info read the
% weblink: http://au.mathworks.com/matlabcentral/answers/17701-rand-state-11
randn('state',SEED1)  %For draws from normal
rand('state',SEED1)   %For draws from uniform

%Create draws
disp('Creating draws.');
% -- Change Start 12/06/2016
% Preparing the Master Draws matrix of dimensions NALTMAX x NP x MASTERDRAWS
MASTERDR=makedraws;   
MASTERDR=permute(MASTERDR,[3,2,1]);   %permute command is used to change the dimensions of a matrix... Read more on http://stackoverflow.com/questions/21100168/how-does-the-permute-function-in-matlab-work
% Draws are generated as NDRAWS X NP X NALTMAX..... 
% Changing the structure of the matrix to NALTMAX X NP X NDRAWS
% DR=reshape(DR,[NALTMAX NP NDRAWS]);
% Repeating the blocks of draws NCSMAX times
% DR=repmat(DR,[NCSMAX 1 1]);
% -- Change End 12/06/2016

% -------- DATA PREPARATION FOR THE MEMBERSHIP MODEL --------
% 
XLCV=zeros(NLCV,NP); % covariates with fixed coefficients.. The size of matrix is NLCV X NP 
%
% The following loop is for the covariates
for n=1:NP;  %loop over people
 xxlcv=XMAT(cp == n, IDLCV(:,1));
 xxlcv=transpose(xxlcv);
 XLCV(:,n)=xxlcv(:,1);
end

% -------- DATA PREPARATION FOR THE CHOICE MODEL --------
%
% Data arrays 
% All variables are differenced from the chosen alternative, i.e. Attribute
% of alternative j - attribute of chosen alternative
% Only nonchosen alternatives are included, since V for chosen alt =0
% This reduces number of calculations for logit prob and eliminates need to
% retain dependent variable.
                                 
XF=zeros(NALTMAX-1,NCSMAX,NF,NP); % Explanatory variables with fixed coefficients 
%                                  for all choice situations, for each person 
S=zeros(NALTMAX-1,NCSMAX,NP); % Identification of the alternatives in each choice situation, for each person

% The following loop is for the fixed variables
for n=1:NP;  %loop over people
 cs=XMAT(cp == n,2);
 yy=XMAT(cp == n,3);
 xxf=XMAT(cp == n, IDF(:,1));
 t1=cs(1,1);
 t2=cs(end,1);
 for t=t1:t2; %loop over choice situations
     k=sum(cs==t)-1; %One less than number of alts = number of nonchosen alts
     S(1:k,1+t-t1,n)=ones(k,1);
        XF(1:k,1+t-t1,:,n)=xxf(cs==t & yy == 0,:)-repmat(xxf(cs==t & yy == 1,:),k,1);
 end
end

% -- Change Start 12/06/2016
% Commneting out this piece of code and inserting the same in prepareX.
% Also changing the clear statements given at last

% X=zeros(NALTMAX-1,NCSMAX,NV,NP,NMEM);  % Defining a coefficient vector
%  
% S=zeros(NALTMAX-1,NCSMAX,NP); % Identification of the alternatives in each choice situation, for each person

% The following loop is for the error components
%for draw=1:NDRAWS;
%for n=1:NP;  %loop over people
% cs=XMAT(cp == n,2);
% yy=XMAT(cp == n,3);
% xx=DR(:,n,draw);
% t1=cs(1,1);
% t2=cs(end,1);
% for t=t1:t2; %loop over choice situations
%     k=sum(cs==t)-1; %One less than number of alts = number of nonchosen alts
%     S(1:k,1+t-t1,n)=ones(k,1);
%     X(1:k,1+t-t1,1,n,draw)=xx(cs==t & yy == 0,:)-repmat(xx(cs==t & yy == 1,:),k,1);
% end
%end
%end

% Calling function prepareX to create random variable matrix from the Master
% Draws matrix. This function will be executed only ONCE... previously, the
% draws were being made in a sequential order which lead to explosion of
% SE. This time I asm selecting draws from the Master Draws matrix at
% random for every individual. If we put this function inside fmincon, the
% LL function will keep on changing at every iteration and thus we won't
% get any meaningful results. Refer to my notes for more information..
X=prepareX;
% -- Change End 12/06/2016

clear global XMAT cp
clear cs yy t1 t2 xx xxf xxlcv k nn


% --------------- NAME THE PARAMETERS TO BE ESTIMATED ---------------------
% Names of parameters for the study are listed in this array
NAMES = {'B_TT' 'B_TTS' 'B_SnGo' 'B_VRC' 'Sigma'};

% Names of covariate parameters for the study. It even includes the class
% specific constant 
NAMESLCV = {'B_CONST' 'B_GENDER' 'B_INCOMEBELOW25' 'B_AGEBELOW40'};
%--------------------------------------------------------------------------

for class = 1:NOC
    for m = 1:NAT
        % Set of all fixed parameters to be estimated
        lc_names_fixed{(class-1)*NAT + m,1} = [NAMES{m} '_' num2str(class)];
    end
    
    for m = 1:NV
        % Set of all random parameters to be estimated
        lc_names_rand{class-1+m,1} = [NAMES{NAT+m} '_' num2str(class)];
    end
    
    if class < NOC  % Keeping the last latent class as refrence
        % Set of all latent class parameters to be estimated
        for m = 1:NLCV
            lc_names{(class-1)*NLCV + m,1} = [NAMESLCV{m} '_' num2str(class)];
        end
    end
end
 
% BOUNDS.... Should be a vector of size 1 X [NOC*(NF+NV)+(NOC-1)*NLCV]. It
% comprises of two components, Lower and Upper bounds
lb=[-100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100];
ub=[100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100];

% NOTE: A good dea to discontinue generating random starting values..
% Prefer giving the values manually for the covariate case
%
%b0 = [rand(SETS,(length(lc_names_fixed) + length(lc_names_rand) + length(lc_names)))];	% This is an array of all parameters for different starting values. This should go into fminunc function. NOTE: 1 row at a time though..
%
% 2 Latent class template
%b0 = [-0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 0.1 0.1 0.1 0.1];     

% 3 Latent class template
b0 = [-0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];

% 4 Latent class template
%b0 = [-0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];

% 5 Latent class template
%b0 = [-0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];

% 6 Latent class template
%b0 = [-0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];

% 7 Latent class template
%b0 = [-0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 -0.1 -0.1 -0.1 -0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];

% b0 and lc_names goes as arguments into estimate.m

disp('Start estimation');
disp('The negative of the log-likelihood is minimized,');
disp('which is the same as maximizing the log-likelihood.');
tic;

for set = 1:SETS
    % Assign one set of starting values for likelihood estimation at a time
    param=b0(set,:);
    %
    options=optimset('LargeScale','off','Display','iter','GradObj','off',...
        'MaxFunEvals',10000,'MaxIter',MAXITERS,'TolX',PARAMTOL,'TolFun',LLTOL,'DerivativeCheck','off');
    [paramhat,fval,exitflag,output,lamda,grad,hessian]=fmincon(@loglik,param,[],[],[],[],lb,ub,@confun,options);    %@ is a handle that passes one function into another function.. Here loglik is a function defined in loglik.m.. Read more on http://au.mathworks.com/help/matlab/matlab_prog/symbol-reference.html
    
    disp(' ');
    disp(['Estimation took ' num2str(toc./60) ' minutes.']);
    disp(' ');
    if exitflag == 1
        disp('Convergence achieved.');
    elseif exitflag == 2
        disp('Convergence achieved by criterion based on change in parameters.');
        if size(PARAMTOL,1)>0
            disp(['Parameters changed less than PARAMTOL= ' num2str(PARAMTOL)]);
        else
            disp('Parameters changed less than PARAMTOL=0.000001, set by default.');
        end
        disp('You might want to check whether this is actually convergence.');
        disp('The gradient vector is');
        grad
    elseif exitflag == 3
        disp('Convergence achieved by criterion based on change in log-likelihood value.');
        if size(PARAMTOL,1)>0
            disp(['Log-likelihood value changed less than LLTOL= ' num2str(LLTOL)]);
        else
            disp('Log-likelihood changed less than LLTOL=0.000001, set by default.');
        end
        disp('You might want to check whether this is actually convergence.');
        disp('The gradient vector is');
        grad
    else
        disp('Convergence not achieved.');
        disp('The current value of the parameters and hessian');
        disp('can be accesses as variables paramhat and hessian.');
        disp('Results are not printed because no convergence.');
        return
    end
    
    disp(['Value of the log-likelihood function at convergence: ' num2str(-fval)]);
    
    % Adding the logic to evaluate AIC and BIC values wrt the final
    % log-likelihood value and the number of parameters
    numParam=length(paramhat);  % Determining the number of estimated parameters
    [aic,bic] = aicbic(-1*fval,numParam,NROWS);   % -fval is the final log-likelihood value.... BIC calculation requires the sample size too... Lower the AIC/BIC value, better is the model
    disp(' ');
    disp('                      AIC                      BIC');
    disp('              ------------------   -----------------------');
    fprintf('                  %10.4f               %10.4f\n', aic,bic);
    disp(' ');
    
    %Calculate standard errors of parameters
    disp(' ');
    disp('Taking inverse of hessian for standard errors.');
    disp(' ');
    ihess=inv(hessian);
    ihess1 = inv(NP*hessian);
    grad1 = NP*grad;
    Jacobian = (grad1*grad1');
    covbhh = ihess*Jacobian*ihess;
    st_final = sqrt(diag(covbhh));
    stderr=sqrt(diag(ihess));
    disp(['The value of grad*inv(hessian)*grad is: ' num2str(grad'*ihess*grad)]);
    
    % Finding the estimated parameters for each latent class (betas, sigma and latent class prevalences)
    for class = 1:NOC
        if NF>0
            if class == 1
                fhat=paramhat(1,(class*NV + (class-1)*NF):(class*NF + (class-1)*NV));   % Change this... paramhat is a row vector
                fsd=stderr((class*NV + (class-1)*NF):(class*NF + (class-1)*NV),1);
            else
                fhat=horzcat(fhat, paramhat(1,(class*NV + (class-1)*NF):(class*NF + (class-1)*NV)));
                fsd=vertcat(fsd, stderr((class*NV + (class-1)*NF):(class*NF + (class-1)*NV),1));
            end
        end
        
        if NV>0
            if class == 1
               bhat=0;  % Setting the mean of EC as zero
               bsd=0;
               what=paramhat(1,class*(NF+NV));
               wsd=stderr(class*(NF+NV),1);
            else
               what=horzcat(what, paramhat(1,class*(NF+NV)));
               wsd=vertcat(wsd, stderr(class*(NF+NV),1));
            end
        end
        
        if class == NOC
            %if class == 1
                lchat=paramhat(1,class*(NF+NV)+1:end);
                lcsd=stderr(class*(NF+NV)+1:end,1);
            %end
            %lchat=vertcat(lchat, paramhat(class*(NF+NV)+1:end,1));
            %lcsd=vertcat(lcsd, stderr(class*(NF+NV)+1:end,1));
        end
    end
    % Taking transpose of standard error vectors
    fsd=transpose(fsd);
    wsd=transpose(wsd);
    lcsd=transpose(lcsd);
    
    disp('RESULTS');
    disp(' ');
    disp(' ')
    if NF>0
        disp('LATENT CLASS COEFFICIENTS');
        disp(' ');
        disp('                      LC      ');
        disp('              ------------------ ');
        disp('                Est         SE ');
        for r=1:length(lc_names);
            fprintf('%-10s %10.4f %10.4f\n', lc_names{r,1}, [lchat(1,r) lcsd(1,r)]);
        end
        disp(' ');
    end
    
    disp(' ');
    disp(' ')
    if NF>0
        disp('FIXED COEFFICIENTS');
        disp(' ');
        disp('                      F      ');
        disp('              ------------------ ');
        disp('                Est         SE ');
        for r=1:length(lc_names_fixed);
            fprintf('%-10s %10.4f %10.4f\n', lc_names_fixed{r,1}, [fhat(1,r) fsd(1,r)]);
        end
        disp(' ');
    end
    
    disp(' ');
    if NV>0;
        disp('RANDOM COEFFICIENTS');
        
        disp(' ');
        disp('                      B                      W');
        disp('              ------------------   -----------------------');
        disp('                 Est     SE            Est         SE');
        for r=1:length(lc_names_rand);
            fprintf('%-10s %10.4f %10.4f %10.4f %10.4f\n', lc_names_rand{r,1}, [bhat bsd what(1,r) wsd(1,r)]);
        end
        
        
        %Create draws of coefficients from B-hat and W-hat
        % NEERAJ SAXENA: I am deleting this block of code as I am no longer
        % doing an RPL... For EC, I don't wish to find the proportion of +ve
        % attribute effects... You can still copy the code from Prof.
        % Train's code from the file "doit.m"
        
        disp(' ');
        
    end
    
    disp(' ');
    disp('ESTIMATED PARAMETERS AND FULL COVARIANCE MATRIX.');
    disp('The estimated values of the parameters are:');
    paramhat
    disp('The covariance matrix for these parameters is:');
    ihess
    
    disp(' ');
    disp('You can access the estimated parameters as variable paramhat,');
    disp('the gradient of the negative of the log-likelihood function as variable grad,');
    disp('the hessian of the negative of the log-likelihood function as variable hessian,');
    disp('and the inverse of the hessian as variable ihess.');
    disp('The hessian is calculated by the BFGS updating procedure.');
end
