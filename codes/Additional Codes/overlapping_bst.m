function [Zb,idx] = overlapping_bst(Z,b)
% PURPOSE: Overlapping Block Bootstrap for a vector time series
% ------------------------------------------------------------
% SYNTAX: Zb = overlappingBB(Z,b);
% ------------------------------------------------------------
% OUTPUT: Zb : (k*b)xkz resampled time series (with k=[n/b])
% ------------------------------------------------------------
% INPUT:  Z  :  nxkz --> vector time series to be resampled
%         b  :  1x1  --> block size (b>=1)
%         If b=1 the Efron's standard iid bootstrap is applied
% ------------------------------------------------------------
% REFERENCES: Knsch, H.R.(1989) "The jacknife and the bootstrap 
% for general stationary observations",The Annals of Statistics, 
% vol. 17, n. 3, p. 1217-1241.
% Davison, A.C. y  Hinkley, D.V. (1997) "Bootstrap methods and 
% their application", Ch. 8: Complex Dependence, Cambridge 
% University Press, Cambridge. U.K.
% ------------------------------------------------------------

% written by:
%  Enrique M. Quilis
%  Macroeconomic Research Department
%  Fiscal Authority for Fiscal Responsibility (AIReF)
%  <enrique.quilis@airef.es>

% Version 1.1 [October 2015]

% adapted by:
% Andrej Sokol
% Research and Strategy Team, GSID, International Directorate
% Bank of England,
% andrej.sokol@bankofengland.co.uk

% ============================================================
% Dimension of time series to be bootstrapped
[n,~] = size(Z);

% Number of blocks
k = fix(n/b);

% Index selection
startI = round(1+(n-b)*rand(1,k));
bstIdx = zeros(k,b);
bstIdx(:,1) = startI;

for bb = 2:b
   
    bstIdx(:,bb) = bstIdx(:,bb-1)+1;
    
end

bstIdx = bstIdx';
idx = bstIdx(:);

Zb = Z(idx,:);

