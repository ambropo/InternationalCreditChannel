function Fcomp = Companion(F,det)
% =======================================================================
% Compute the companion matrix for a VAR estimated with VARmodel. 
% =======================================================================
% Fcomp = Companion(F,det)
% -----------------------------------------------------------------------
% INPUT
%   - F: matrix of coefficients from VARmodel
%   - det: number of deterministic components (see VARmodel)
% -----------------------------------------------------------------------
% OUTPUT
%   - Fcomp: companion matrix
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com


%% Check inputs
%===============================================
if ~exist('det','var')
    det = 0;
end
if size(F,2)/size(F,1)<1
    F = F';
end

%% Retrieve parameters and preallocate variables
%===============================================
F = F(:,1+det:end);
nvar = size(F,1);
nlags = size(F,2)/size(F,1);
Fcomp = [F(:,1:nvar*nlags) ; eye(nvar*(nlags-1)) zeros(nvar*(nlags-1),nvar)];
