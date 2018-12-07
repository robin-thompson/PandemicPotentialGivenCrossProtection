clc
clear all
close all

%% This code accompanies the manuscript entitled "Increased frequency of travel in the
%% presence of cross-immunity may act to decrease the chance of a global pandemic"
%% by Thompson et al. For further information about the paper or this code, please email
%% robin.thompson@chch.ox.ac.uk

%% We request that users cite the original publication when referring to this code or
%% any results generated from it.

%% Please note - this code may take some time to run.

% Here we estimate three quantities using analytic calculations (along with numerical solution of ODEs) for the full two population model:
% 1. The probability of a major epidemic of the HV strain
% 2. The mean final size of major epidemics of the HV strain
% 3. The mean final size of all outbreaks of the HV strain (minor outbreaks and major epidemics)

% These results correspond to e.g. individual tiles in the left column of
% Fig S1 of our paper.


global beta
global alpha
global mu

% Set the value of the between-patch travel parameter, lambda
lambda = 2.5*10^(-4);
% Set the value of the within-patch travel parameter, kappa
kappa = 1.5;


% Set up other parameter values
N = 1000;
muNP = 1/7;
muP = 1/7;
RzNP = kappa*2;
RzP = kappa*3;
alpha = 0.7;

probVector = zeros(8,1); %Vector for probability of each of cases 1,2,3,...8 occurring (for descriptions of each case, please see our manuscript)
finalSizesVector = zeros(8,1);

% Calculate final size of an LV strain major epidemic in a single patch
% Note, below the notation "NP" refers to the
% non-pandemic - as opposed to the pandemic, or "P" - strain.
% LV strain = NP strain
% HV strain = P strain
Sz = N;
syms x
eqn = x == N - Sz*exp(-x*RzNP/N);
solx = solve(eqn,x);
RinfNP = eval(solx);

% Calculate final size of an HV strain major epidemic in a single patch, if
% the LV strain has not previously invaded that patch
Sz = N;
syms x
eqn = x == N - Sz*exp(-x*RzP/N);
solx = solve(eqn,x);
RinfP = eval(solx);

% Calculate effective reproduction number of the HV strain in a single patch, if
% the LV strain has previously invaded that patch
ReffP = RzP*(1 - alpha*RinfNP/N);

% Numerically find final size of an HV strain major epidemic in a single patch, if
% the LV strain has previously invaded that patch
mu = muP;
beta = RzP*mu/N;
fracProtected = RinfNP/N;
tVals = [0:0.01:10000];
I = 1;
R = 0;
S = N - I - R;
SI = floor(fracProtected*S); % S immune
SN = S - SI; % S not immune
yo = [SI SN I R];
[t y] = ode45('modelDEsFunction',tVals, yo);

RinfeffP = y(length(y(:,1)),4);

% Calculate the probabilities of each of the cases described in the
% manuscript, and their corresponding final sizes
if ReffP > 1
    
    probVector(1) = (1-(lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP)*(1/ReffP);
    probVector(2) = (1-(lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP)*0.5*(1 - 1/ReffP)*((lambda/(ReffP*(lambda + muP)) + muP/(lambda + muP))^RinfeffP);
    probVector(3) = probVector(2);
    probVector(4) = (1-(lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP)*(1 - 1/ReffP)*(1 - (lambda/(ReffP*(lambda + muP)) + muP/(lambda + muP))^RinfeffP);
    probVector(5) = ((lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP)*(0.5*(1/ReffP) + 0.5*(1/RzP));
    probVector(6) = ((lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP)*0.5*(1 - 1/ReffP)*((lambda/(RzP*(lambda + muP)) + muP/(lambda + muP))^RinfeffP);
    probVector(7) = ((lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP)*0.5*(1 - 1/RzP)*((lambda/(ReffP*(lambda + muP)) + muP/(lambda + muP))^RinfP);
    probVector(8) = ((lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP)*(0.5*(1 - 1/ReffP)*(1 - (lambda/(RzP*(lambda + muP)) + muP/(lambda + muP))^RinfeffP) + 0.5*(1 - 1/RzP)*((1 - (lambda/(ReffP*(lambda + muP)) + muP/(lambda + muP))^RinfP)));
    
else
    
    probVector(1) = (1-(lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP);
    probVector(2) = 0;
    probVector(3) = 0;
    probVector(4) = 0;
    probVector(5) = ((lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP)*(0.5 + 0.5*(1/RzP));
    probVector(6) = 0;
    probVector(7) = ((lambda/(RzNP*(lambda + muNP)) + muNP/(lambda + muNP))^RinfNP)*(0.5*(1-1/RzP))
    probVector(8) = 0;
    
end

finalSizesVector(1) = 0;
finalSizesVector(2) = RinfeffP;
finalSizesVector(3) = RinfeffP;
finalSizesVector(4) = 2*RinfeffP;
finalSizesVector(5) = 0;
finalSizesVector(6) = RinfeffP;
finalSizesVector(7) = RinfP;
finalSizesVector(8) = RinfeffP + RinfP;

sumProbs = sum(probVector);
probMO = probVector(2) + probVector(3) + probVector(4) + probVector(6) + probVector(7) + probVector(8);
expectedFinalSizeGivenMO = (probVector(2).*finalSizesVector(2) + probVector(3).*finalSizesVector(3) + probVector(4).*finalSizesVector(4) + probVector(6).*finalSizesVector(6) + probVector(7).*finalSizesVector(7) + probVector(8).*finalSizesVector(8))/(probVector(2) + probVector(3) + probVector(4) + probVector(6) + probVector(7) + probVector(8));
expectedFinalSize = sum(probVector.*finalSizesVector);

clc

probMO
expectedFinalSizeGivenMO
expectedFinalSize