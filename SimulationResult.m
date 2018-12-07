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

% Here we estimate three quantities using stochastic simulations of the full two population model:
% 1. The probability of a major epidemic of the HV strain
% 2. The mean final size of major epidemics of the HV strain
% 3. The mean final size of all outbreaks of the HV strain (minor outbreaks and major epidemics)

% These results correspond to e.g. individual tiles in the right column of
% Fig S1 of our paper.



% Set the value of the between-patch travel parameter, lambda
lambda = 2.5*10^(-4);
% Set the value of the within-patch travel parameter, kappa
kappa = 1.5;


% Set up other parameter values
N = 1000;
mu = 1/7;
RzNP = kappa*2;
RzP = kappa*3;
alpha = 0.7;
MOThreshold = 20;

probVector = zeros(8,1); %Vector for probability of each of cases 1,2,3,...8 occurring (for descriptions of each case, please see our manuscript)
finalSizesVector = zeros(8,1);
nTimesVector = zeros(8,1);

nSims = 10000;

simNo = 1;
while simNo < nSims + 0.5
    
    % First, run the LV strain outbreak (note, below the notation "NP" refers to the
    % non-pandemic - as opposed to the pandemic, or "P" - strain)
    % LV strain = NP strain
    % HV strain = P strain
    
    SNP = zeros(2,1);
    INP = zeros(2,1);
    RNP = zeros(2,1);
    
    SNP(1) = N - 1;
    INP(1) = 1;
    RNP(1) = 0;
    SNP(2) = N;
    INP(2) = 0;
    RNP(2) = 0;
    
    beta = RzNP*mu/N;
    
    while sum(INP) > 0
        
        alphaZero = beta*INP(1)*SNP(1) + beta*INP(2)*SNP(2) + mu*INP(1) + mu*INP(2) + lambda*(sum(SNP) + sum(INP) + sum(RNP));
        
        r = rand();
        if r < beta*INP(1)*SNP(1)/alphaZero
            eventType = 1;
        end
        if ((r > beta*INP(1)*SNP(1)/alphaZero) && (r < (beta*INP(1)*SNP(1) + beta*INP(2)*SNP(2))/alphaZero))
            eventType = 2;
        end
        if ((r > (beta*INP(1)*SNP(1) + beta*INP(2)*SNP(2))/alphaZero) && (r < (beta*INP(1)*SNP(1) + beta*INP(2)*SNP(2) + mu*INP(1))/alphaZero))
            eventType = 3;
        end
        if ((r < (beta*INP(1)*SNP(1) + beta*INP(2)*SNP(2) + mu*INP(1) + mu*INP(2))/alphaZero) && (r > (beta*INP(1)*SNP(1) + beta*INP(2)*SNP(2) + mu*INP(1))/alphaZero))
            eventType = 4;
        end
        if ((r > (beta*INP(1)*SNP(1) + beta*INP(2)*SNP(2) + mu*INP(1) + mu*INP(2))/alphaZero) && (r < (beta*INP(1)*SNP(1) + beta*INP(2)*SNP(2) + mu*INP(1) + mu*INP(2) + lambda*(sum(SNP) + sum(INP) + sum(RNP)))/alphaZero))
            eventType = 5;
        end
        
        if eventType == 1
            SNP(1) = SNP(1) - 1;
            INP(1) = INP(1) + 1;
        end
        if eventType == 2
            SNP(2) = SNP(2) - 1;
            INP(2) = INP(2) + 1;
        end
        if eventType == 3
            INP(1) = INP(1) - 1;
            RNP(1) = RNP(1) + 1;
        end
        if eventType == 4
            INP(2) = INP(2) - 1;
            RNP(2) = RNP(2) + 1;
        end
        if eventType == 5
            rWO = rand();
            if rWO < (sum(SNP)/sum(SNP + INP + RNP))
                whichMove = 1;
            end
            if ((rWO > (sum(SNP)/sum(SNP + INP + RNP))) && (rWO < ((sum(SNP) + sum(INP))/sum(SNP + INP + RNP))))
                whichMove = 2;
            end
            if  (rWO > ((sum(SNP) + sum(INP))/sum(SNP + INP + RNP)))
                whichMove = 3;
            end
            
            if whichMove == 1
                rWHICH = rand();
                if rWHICH < SNP(1)/sum(SNP)
                    SNP(1) = SNP(1) - 1;
                    SNP(2) = SNP(2) + 1;
                else
                    SNP(1) = SNP(1) + 1;
                    SNP(2) = SNP(2) - 1;
                end
            end
            if whichMove == 2
                rWHICH = rand();
                if rWHICH < INP(1)/sum(INP)
                    INP(1) = INP(1) - 1;
                    INP(2) = INP(2) + 1;
                else
                    INP(1) = INP(1) + 1;
                    INP(2) = INP(2) - 1;
                end
            end
            if whichMove == 3
                rWHICH = rand();
                if rWHICH < RNP(1)/sum(RNP)
                    RNP(1) = RNP(1) - 1;
                    RNP(2) = RNP(2) + 1;
                else
                    RNP(1) = RNP(1) + 1;
                    RNP(2) = RNP(2) - 1;
                end
            end
            
        end
        
    end
    
    % ensure that the LV strain outbreak was a major epidemic
    if RNP(1) > MOThreshold
        
        % Run HV strain epidemic
        SPNonImmune = zeros(2,1);
        SPImmune = zeros(2,1);
        IP = zeros(2,1);
        RP = zeros(2,1);
        
        r = rand();
        if r < 0.5
            IP(1) = 1;
        else
            IP(2) = 1;
        end
        SP1 = N - RP(1) - IP(1);
        SP2 = N - RP(2) - IP(2);
        
        SPImmune(1) = round(SP1*(RNP(1)/N));
        SPNonImmune(1)= SP1 - SPImmune(1);
        
        SPImmune(2) = round(SP2*(RNP(2)/N));
        SPNonImmune(2)= SP2 - SPImmune(2);
        
        beta = RzP*mu/N;
        
        while sum(IP) > 0
            
            eventRates = zeros(14,1);
            eventRates(1) = beta*(1 - alpha)*IP(1)*SPImmune(1);
            eventRates(2) = beta*IP(1)*SPNonImmune(1);
            eventRates(3) = mu*IP(1);
            eventRates(4) = beta*(1 - alpha)*IP(2)*SPImmune(2);
            eventRates(5) = beta*IP(2)*SPNonImmune(2);
            eventRates(6) = mu*IP(2);
            eventRates(7) = lambda*SPNonImmune(1);
            eventRates(8) = lambda*SPImmune(1);
            eventRates(9) = lambda*IP(1);
            eventRates(10) = lambda*RP(1);
            eventRates(11) = lambda*SPNonImmune(2);
            eventRates(12) = lambda*SPImmune(2);
            eventRates(13) = lambda*IP(2);
            eventRates(14) = lambda*RP(2);
            
            alphaZero = sum(eventRates);
            
            r = rand();
            if r < eventRates(1)/alphaZero
                eventType = 1;
            end
            if ((r > eventRates(1)/alphaZero) && (r < sum(eventRates(1:2))/alphaZero))
                eventType = 2;
            end
            if ((r > sum(eventRates(1:2))/alphaZero) && (r < sum(eventRates(1:3))/alphaZero))
                eventType = 3;
            end
            if ((r > sum(eventRates(1:3))/alphaZero) && (r < sum(eventRates(1:4))/alphaZero))
                eventType = 4;
            end
            if ((r > sum(eventRates(1:4))/alphaZero) && (r < sum(eventRates(1:5))/alphaZero))
                eventType = 5;
            end
            if ((r > sum(eventRates(1:5))/alphaZero) && (r < sum(eventRates(1:6))/alphaZero))
                eventType = 6;
            end
            if ((r > sum(eventRates(1:6))/alphaZero) && (r < sum(eventRates(1:7))/alphaZero))
                eventType = 7;
            end
            if ((r > sum(eventRates(1:7))/alphaZero) && (r < sum(eventRates(1:8))/alphaZero))
                eventType = 8;
            end
            if ((r > sum(eventRates(1:8))/alphaZero) && (r < sum(eventRates(1:9))/alphaZero))
                eventType = 9;
            end
            if ((r > sum(eventRates(1:9))/alphaZero) && (r < sum(eventRates(1:10))/alphaZero))
                eventType = 10;
            end
            if ((r > sum(eventRates(1:10))/alphaZero) && (r < sum(eventRates(1:11))/alphaZero))
                eventType = 11;
            end
            if ((r > sum(eventRates(1:11))/alphaZero) && (r < sum(eventRates(1:12))/alphaZero))
                eventType = 12;
            end
            if ((r > sum(eventRates(1:12))/alphaZero) && (r < sum(eventRates(1:13))/alphaZero))
                eventType = 13;
            end
            if (r > sum(eventRates(1:13))/alphaZero)
                eventType = 14;
            end
            
            if eventType == 1
                SPImmune(1) = SPImmune(1) - 1;
                IP(1) = IP(1) + 1;
            end
            if eventType == 2
                SPNonImmune(1) = SPNonImmune(1) - 1;
                IP(1) = IP(1) + 1;
            end
            if eventType == 3
                IP(1) = IP(1) - 1;
                RP(1) = RP(1) + 1;
            end
            if eventType == 4
                SPImmune(2) = SPImmune(2) - 1;
                IP(2) = IP(2) + 1;
            end
            if eventType == 5
                SPNonImmune(2) = SPNonImmune(2) - 1;
                IP(2) = IP(2) + 1;
            end
            if eventType == 6
                IP(2) = IP(2) - 1;
                RP(2) = RP(2) + 1;
            end
            if eventType == 7
                SPNonImmune(1) =  SPNonImmune(1) - 1;
                SPNonImmune(2) =  SPNonImmune(2) + 1;
            end
            if eventType == 8
                SPImmune(1) =  SPImmune(1) - 1;
                SPImmune(2) =  SPImmune(2) + 1;
            end
            if eventType == 9
                IP(1) = IP(1) - 1;
                IP(2) = IP(2) + 1;
            end
            if eventType == 10
                RP(1) = RP(1) - 1;
                RP(2) = RP(2) + 1;
            end
            if eventType == 11
                SPNonImmune(2) =  SPNonImmune(2) - 1;
                SPNonImmune(1) =  SPNonImmune(1) + 1;
            end
            if eventType == 12
                SPImmune(2) =  SPImmune(2) - 1;
                SPImmune(1) =  SPImmune(1) + 1;
            end
            if eventType == 13
                IP(2) = IP(2) - 1;
                IP(1) = IP(1) + 1;
            end
            if eventType == 14
                RP(2) = RP(2) - 1;
                RP(1) = RP(1) + 1;
            end
            
        end
        
        
        if ((RNP(2) >= MOThreshold) && ((RP(1) < MOThreshold) && (RP(2) < MOThreshold)))
            nTimesVector(1) = nTimesVector(1) + 1;
            finalSizesVector(1) = finalSizesVector(1) + RP(1) + RP(2);
        end
        if ((RNP(2) >= MOThreshold) && ((RP(1) >= MOThreshold) && (RP(2) < MOThreshold)))
            nTimesVector(2) = nTimesVector(2) + 1;
            finalSizesVector(2) = finalSizesVector(2) + RP(1) + RP(2);
        end
        if ((RNP(2) >= MOThreshold) && ((RP(1) < MOThreshold) && (RP(2) >= MOThreshold)))
            nTimesVector(3) = nTimesVector(3) + 1;
            finalSizesVector(3) = finalSizesVector(3) + RP(1) + RP(2);
        end
        if ((RNP(2) >= MOThreshold) && ((RP(1) >= MOThreshold) && (RP(2) >= MOThreshold)))
            nTimesVector(4) = nTimesVector(4) + 1;
            finalSizesVector(4) = finalSizesVector(4) + RP(1) + RP(2);
        end
        if ((RNP(2) < MOThreshold) && ((RP(1) < MOThreshold) && (RP(2) < MOThreshold)))
            nTimesVector(5) = nTimesVector(5) + 1;
            finalSizesVector(5) = finalSizesVector(5) + RP(1) + RP(2);
        end
        if ((RNP(2) < MOThreshold) && ((RP(1) >= MOThreshold) && (RP(2) < MOThreshold)))
            nTimesVector(6) = nTimesVector(6) + 1;
            finalSizesVector(6) = finalSizesVector(6) + RP(1) + RP(2);
        end
        if ((RNP(2) < MOThreshold) && ((RP(1) < MOThreshold) && (RP(2) >= MOThreshold)))
            nTimesVector(7) = nTimesVector(7) + 1;
            finalSizesVector(7) = finalSizesVector(7) + RP(1) + RP(2);
        end
        if ((RNP(2) < MOThreshold) && ((RP(1) >= MOThreshold) && (RP(2) >= MOThreshold)))
            nTimesVector(8) = nTimesVector(8) + 1;
            finalSizesVector(8) = finalSizesVector(8) + RP(1) + RP(2);
        end
        
        simNo = simNo + 1;
        
    end
    
end


% Use the simulation results to calculate the quantities of interest
probVector = nTimesVector./nSims; 

for i = 1:length(finalSizesVector)
    if nTimesVector(i) > 0
        finalSizesVector(i) = finalSizesVector(i)/nTimesVector(i);
    end
end

sumProbs = sum(probVector);
probMO = probVector(2) + probVector(3) + probVector(4) + probVector(6) + probVector(7) + probVector(8);
expectedFinalSizeGivenMO = (probVector(2).*finalSizesVector(2) + probVector(3).*finalSizesVector(3) + probVector(4).*finalSizesVector(4) + probVector(6).*finalSizesVector(6) + probVector(7).*finalSizesVector(7) + probVector(8).*finalSizesVector(8))/(probVector(2) + probVector(3) + probVector(4) + probVector(6) + probVector(7) + probVector(8));
expectedFinalSize = sum(probVector.*finalSizesVector);

clc

probMO
expectedFinalSizeGivenMO
expectedFinalSize