% Financial shocks, credit spreads, and the international credit channel
% Ambrogio Cesa-Bianchi and Andrej Sokol, July 2021
%
% Replication file for all VAR-based figures in the paper, appendix and
% supplement.


%% Housekeeping
clear all; clear session; close all; clc
warning off all
fontSize   = 13;
highQuality = 0;
folderOut = '..\figures\';
frequency = 'm';
addpath(genpath('..\codes'));
addpath('..\scripts');

% Initialize
userFO = '1979m7';
userLO = '2016m12';
sNames = {'Monetary Policy';'CB Information';'Financial';'Demand';'Supply'};
dynamics = 1; % set to zero to get no influence of US dynamics on UK / set to 1 for baseline
mkdir(folderOut)
fname = 'VAR_data_2019m12.xlsx';

% VAR specification
nLag = 12; nLagSOE = 12; nLagXSOE = 3; VARCase = 1; ident='adhoc'; VARtol=0.999; NOBStol=30; nSteps=36; pctg=90; pctg2 = 68; nDraws=2e6;
nShocks = 5; nPeriods = 3; maxIters = 1e4; modUnc = 1; finShockSignOnly = 1;
pctg_inf = (100-pctg)/2;
pctg_sup = 100 - (100-pctg)/2;
pctg_inf2 = (100-pctg2)/2;
pctg_sup2 = 100 - (100-pctg2)/2;


%% US VAR (BASELINE)
STEP1_US_VAR_baseline;


%% UK VAR (BASELINE AND COUNTERFACTUAL)
STEP2_UK_VAR_BASELINE_AND_CF;


%% Baseline Plots (Figures 4-9 and A1a)
STEP3_BASELINE_PLOTS;


%% Counterfactual Plots (Figures 10 and B4-B9)
STEP4_COUNTERFACTUAL_PLOTS;


%% UK Inflation Targeting Sample (Figure S2)
STEP5_UK_IT_SAMPLE;


%% UK VAR (OTHER SPREADS, Figure 11)
STEP6_UK_VAR_MORTGAGE_SPR;
STEP7_UK_VAR_LIBOR_SPR;


%% Long rates spread (Figure A1b)
STEP8_10y_SPR;


%% Additional shock (Figure S1)
STEP9_ADDITIONAL_SHOCK;


%% Pre-2008 Sample (Figure S3)
STEP10_PRE_GFC;


%% Unemployment (Figure S4)
STEP11_UNEMPLOYMENT;


