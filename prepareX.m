% ----------------------------------------------------------------
% 12/06/2016: Newly created function prepareX which will prepare the random
% variable matrix X for further computation. The code is quite similar to
% the one in doitr previously. IDEA is to reduce the noise effect by randomly
% selecting a draw matrix from Master Draws. The function will be called
% only ONCE during the simulation.
% ----------------------------------------------------------------

function x=prepareX
 
global NALTMAX NCSMAX NP NV NMEM NDRAWS MASTERDR AMPL MASTERDRAWS XMAT cp

% Define an empty array to populate the random variable
x=zeros(NALTMAX-1,NCSMAX,NV,NP,NMEM);  
% Define an empty smaller matrix to store draws from the Master Draws
DR=zeros(NALTMAX,NP,NDRAWS);

% FIRST step is to extract a smaller matrix from Master Draws 
% Size of master Draws  -- NALTMAXxNPxMASTERDRAWS
% Size of smaller Draws -- NALTMAXxNPxNDRAWS
% LOGIC for extraction: Pick Ndraws values at random from the Master Draws for every individual

for n=1:NP 
    % Select a value at random within the range [1,(AMPL-1)*NDRAWS]. This
    % marks the starting index for the NDRAWS to be selected. The Endig
    % index will be equal to the starting index + (NDRAWS-1). This way we
    % selected NDRAWS for an individual n
    if AMPL > 1
       startind=ceil(rand(1)*(MASTERDRAWS - NDRAWS - 1));
    else
       startind=1; 
    end
    endind=startind+(NDRAWS-1);
    DR(:,n,:)=MASTERDR(:,n,startind:endind);
end

% SECOND step is to populate random variable matrix using the smaller DR
% matrix. For doing that, the previously existing code will be used. So I
% commented out the relevant code from doitr.m and pasted the same in this
% function.
% Following is the code copied from doitr.m

% Repeating the blocks of draws NCSMAX times. The dimension of DR now
% becomes (NCSMAX X NALTMAX) x NP x NDRAWS
DR=repmat(DR,[NCSMAX 1 1]);

S=zeros(NALTMAX-1,NCSMAX,NP); % Identification of the alternatives in each choice situation, for each person

% The following loop is for the error components
for draw=1:NDRAWS;
    for n=1:NP;  %loop over people
        cs=XMAT(cp == n,2);
        yy=XMAT(cp == n,3);
        xx=DR(:,n,draw);
        t1=cs(1,1);
        t2=cs(end,1);
        for t=t1:t2; %loop over choice situations
            k=sum(cs==t)-1; %One less than number of alts = number of nonchosen alts
            S(1:k,1+t-t1,n)=ones(k,1);
            x(1:k,1+t-t1,1,n,draw)=xx(cs==t & yy == 0,:)-repmat(xx(cs==t & yy == 1,:),k,1);
        end
    end
end
