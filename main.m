% -----------------------------------------------------------------
%     LATENT CLASS ANALYSIS - CLASS SPECIFICATION WITH COVARIATES
%           Written by: Neeraj Saxena, 20th May 2016 
%            emailid: n.saxena@student.unsw.edu.au
%                 !! SHRI GANESHAAYE NAMAH !!
% -----------------------------------------------------------------
%
% Introduction: This is the main file which starts the LCA process
% This file prepares the data, creates the list of parameters to be
% estimated and set their names. The parameters will be different for each
% class.
% This LCA model considers Error Compnent Logit at the lower level while
% the upper level is expressed in terms of covariates... Each latent class has its own set of beta
% and sigma values in the lower level. For the upper level, the last latent
% class has betas set to ZERO.

% -----------------------------------------------------------------
%                            CHANGE LOG

% 12/06/2016: There seems to be an issue with the way draws are being
% generated at the moment. Make draws is generating draws in a sequential 
% manner. This is probably causing Standard Errors to explode. As a
% possible fix, I'll create a Master Draws matrix, random draws from which
% will be made for every individual. This would ZERO out the the noise
% effect.
%  -----------------------------------------------------------------

clear all

% Declare GLOBAL variables
% GLOBAL variables are all in caps
% DO NOT CHANGE ANY OF THESE 'global' STATEMENTS
global NP NCS NROWS NAT 
global NV NAMES 
global IDF NF 
global IDLCV NLCV
global DRAWTYPE NDRAWS SEED1 
global NALTMAX NCSMAX
global X XF XLCV S DR 
global XMAT
global NMEM NTAKES 
global NOC SETS
% -- Change Start 12/06/2016
global AMPL MASTERDRAWS MASTERDR
global cp
% -- Change End 12/06/2016

% OUTPUT FILE
% Put the name you want for your output file (including full path if not the current 
% working directory) after words "delete" and "diary".
% The 'diary off' and 'delete filename' commands close and delete the previous version 
% of the file created during your current matlab session or in any previous sessions. 
% If you want to append the new output to the old output, then 
% put % in front of the 'diary off' and 'delete filename' commands (or erase them).

diary off
delete myrun249.out
diary myrun249.out

% TITLE
% Put a title for the run in the quotes below, to be printed at the top of the output file.
disp 'Latent Class Model with Covariates & Error Component Logit as the Choice Model.'

% DATA

% Number of people (decision-makers) in dataset.. CHANGE THIS FIELD !!!!
NP=249;        

% Number of choice situations in dataset. This is the number faced by all
% the people combined.. CHANGE THIS FIELD !!!!
NCS=2490;  

% Number of attributes defining an alternative. These attributes will be
% having fixed parameters
NAT=4;

%----------- SET THE COVARIATES ----------
% Number of covariates to be estimated in the lower level. NLCV is only used to define the variable matrix. 
% NOTE: Column 8 is always 1... I am keeping this column just to define
% CLASS SPECIFIC CONSTANTS... 
IDLCV=[8;9;12;19];
NLCV=size(IDLCV,1);

%----------- SET THE ATTRIBUTES FOR CHOICE MODEL ----------
% Number of attributes to be estimated in the lower level. NF is only used to define the variable matrix. 
IDF=[4;5;6;7];
NF=size(IDF,1);

% The error component parameter (sigma). It follows standard ND~(0,1). Also the
% sigma value is same for all alternatives in a class to maintain
% homoscedasticity. Think as if NV is the sigma to be estimated within a single latent class 
NV=1;

% Load and/or create XMAT, a matrix that contains the data.
%
% XMAT must contain one row of data for each alternative in each choice situation for each person.
% The rows are grouped by person, and by choice situations faced by each person.
% The number of rows in XMAT must be NROWS, specified above.
% The columns in XMAT are variable that describe the alternative.
% 
% The *first* column of XMAT identifies the person who faced this alternative. 
% The people must be numbered sequentially from 1 to NP, in ascending order.
% All alternatives for a given person must be grouped together.
% The *second* column of XMAT identifies the choice situation. The choice
% situations must be numbered sequentially from 1 to NCS.
% All alternatives for a given choice situation must be grouped together.
% The *third* column of XMAT identifies the chosen alternatives (1 for
% chosen, 0 for not). One and only one alternative must be chosen for each
% choice situation.
% The remaining columns of XMAT can be any variables.

XMAT=load('SnGo_249 responses.txt');  %The variables are described below

% To help you keep up with the variables, list the variables in XMAT here.
% Start each line with % so that matlab sees that it is a comment rather than a command.
% NOTES for XMAT for Stop-&-go run:
% This dataset is for people's choice among three routes in stated-preference
% experiments. Each person faced with 10 experiments. Each
% experiment contained 3 alternatives (Status Quo + 2 hypothetical) representing three different routes, 
% each defined in terms of 4 attributes (TT, TTS, SnGo and VRC). The person stated which
% of the three routes he/she would take for their travel.
% The variables in XMAT are:
% 1. Person number (1-NP)            MUST BE THIS. DO NOT CHANGE.
% 2. Choice situation number (1-NCS) MUST BE THIS. DO NOT CHANGE.
% 3. Chosen alternative (1/0)        MUST BE THIS. DO NOT CHANGE.
% 4. Total travel time on a route (minutes)
% 5. Time spent in stop-&-go traffic (minutes)
% 6. Number of stop-&-go experienced (number)
% 7. Vehicle running cost for the trip (AU $)
% 8. Alternate specific constant variable which is set to ONE
% 9 onwards is the socio-demographic information of the individual

% Type of draws to use in simulation
% 1=pseudo-random draws
% 2=standard Halton draws
% 3=shifted and shuffled Halton draws
% 4=Modified Latin Hypercube Sampling, shifted and shuffled 
% 5=create your own draws or load draws from file
DRAWTYPE=2;
disp(' ');
switch DRAWTYPE
    case 2
        disp('Random Draws to be made using the Standard Halton process.')
    case 3
        disp('Random Draws to be made using Shuffled and Shifted Halton process.')
    case 4
        disp('Random Draws to be made using MLHS process.')
    otherwise
        disp('Invalid drawing technique.')
end  

% Number of draws to be generated per person in simulation.
NDRAWS=1000;
disp(' ');
fprintf('Number of Random Draws to be made --> %d\n', NDRAWS);

% Memory use
% Give the number of draws that you want held in memory at one time.
% This number must be evenly divisible into the number of draws.
% That is NDRAWS./NMEM must be a positive integer.
% To hold all draws in memory at once, set NMEM=NDRAWS.
% A larger value of NMEM requires fewer reads from disc but 
% uses more memory which can slow-down the calculations and increases 
% the chance of running out of memory.
% If DRAWTYPE=5, then you must set NMEM=NDRAWS
NMEM=NDRAWS;

NTAKES=NDRAWS./NMEM; %This tells code how many passes through the draws are needed 
                           % given NMEM people in each pass.
                           
% -- Change Start 12/06/2016
% The variable AMPL amplifies the size of the Master Draws matrix to be generated. The
% Master Draws matrix size will be [MASTERDRAWS x NP x NALTMAX]
AMPL=10;
disp(' ');
fprintf('Amplification Factor for Master Draw Matrix --> %d\n', AMPL);
MASTERDRAWS=AMPL*NDRAWS;
% -- Change End 12/06/2016


%--------------------------- IMPORTANT NOTE -------------------------------
% This model isn't a well behaved one. The starting value do matter
% critically now. One might see an Inf Hessian matrix upon estimation... 
% REMEDIES:
%     1. Set the seed value carefully
%     2. Add one variable ata time and gradually increment the starting
%     values
%     3. If nothing works, pick up few good starting values and copy them
%     for every other latent class
%
% Refer to readme and my notes for more information
% Set seed for the random number generator.
SEED1 = 28593;
% -------------------------------------------------------------------------

%----------- SET THE NUMBER OF LATENT CLASSSES ----------
% Defining a new variable to set the number of latent classes
NOC=3;
disp(' ');
fprintf('Number of Latent Classes being estimated --> %d\n\n', NOC);

% LC models get often stuck in local minima. Therefore, it is necessary to
% use a series of different starting values 
% Number of starting value sets
SETS = 1;
   
% OPTIMIZATION 
% Maximum number of iterations for the optimization routine.
% The code will abort after ITERMAX iterations, even if convergence has
% not been achieved. The default is 400, which is used when MAXITERS=[];
MAXITERS=[];

% Convergence criterion based on the maximum change in parameters that is considered
% to represent convergence. If all the parameters change by less than PARAMTOL 
% from one iteration to the next, then the code considers convergence to have been
% achieved. The default is 0.000001, which is used when PARAMTOL=[];
PARAMTOL=0.000001;

% Convergence criterion based on change in the log-likelihood that is
% considered to represent convergence. If the log-likelihood value changes
% less than LLTOL from one iteration to the next, then the optimization routine
% considers convergence to have been achieved. The default is 0.000001,
% which is used when LLTOL=[];
LLTOL=[];

%Do not change the next line. It runs the model.
doitr

