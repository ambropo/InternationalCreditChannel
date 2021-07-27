function [PC, TABLE] = PairCorr(DATA,labels)
% =======================================================================
% Compute pairwise correlation of a panel of time series DATA with T 
% observations and N variables. Computes the pairwise correlation of both 
% levels and log-differences of the original series in DATA.
% 
% Note: If a series (ie, column) has a NaN, the whole series is treated as 
% NaN. See PairCorrUnbalanced for pairwise correlations with NaNs
% =======================================================================
% [PC, TABLE] = PairCorr(DATA,labels)
% -----------------------------------------------------------------------
% INPUT
%	- DATA: matrix DATA T (observations) x N (variables)
% -----------------------------------------------------------------------
% OPTIONAL INPUT
%   - labels: names of each variable j. Default "Variable"
% -----------------------------------------------------------------------
% TABLEPUT
%	- PC: matrix of pairwise correlation N (variables) x 2 (levels 
%       and log-diff)
%	- TABLE: formatted table of pairwise correlation with titles
% =========================================================================
% Example
% DATA = rand(50,4);
% [PC TABLE] = PairCorr(DATA)
% =========================================================================
% Ambrogio Cesa Bianchi, March 2015



%% Preliminaries: define the dimension of the matrix of interest
% =========================================================================
[~, col] = size(DATA);

% If no names are provided set it to 'Variable'
if ~exist('labels','var')
    labels(1,1:col) = {'Variable'};
end

% If labels are entered as jx1 vector, transpose it
if size(labels,1) > 1
    labels = labels';
end

% If no option for excel TABLEput is provided, set it to 0 (no TABLEput)
if ~exist('write','var')
    write = 0;
end

%% Compute pairwise correlation
% =========================================================================
% Define the matrix for the levels 
X = DATA;
% If there are columns with NaNs make the whole column NaN
X(:,isnan(X(1,:))==1) = NaN;

% Compute the first differences
x = log(X(2:end,:)) - log(X(1:end-1,:));

%Compute the correlation matrix
X_corr = corrcoef(X);
x_corr = corrcoef(x);

% Find the NaNs in the correlation matrix
nans = isnan(X_corr);

% Define the number of non-NaN elements in the correlation matrix
n = sum(~nans);
n(n==0) = NaN;

% Set the NaNs of the correlation matrix to zero
X_corr(nans) = 0;
x_corr(nans) = 0;

% Compute the pairwise correlation
X_PairCorr = (sum(X_corr) - 1) ./ (n-1);
x_PairCorr = (sum(x_corr) - 1) ./ (n-1);

PC = [X_PairCorr ; x_PairCorr]';

% Table with titles
title = {'' , 'Level', 'First Diff.'};
TABLE = [labels  ; num2cell(X_PairCorr) ; num2cell(x_PairCorr)];
TABLE = [title ; TABLE'];