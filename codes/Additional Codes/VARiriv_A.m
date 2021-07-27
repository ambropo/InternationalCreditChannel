function [IRF, IVout] = VARiriv_A(VAR,VARopt,IV,varargin)
% =======================================================================
% Compute IRFs for a VAR model estimated with VARmodel using the external 
% instrument approach of Stock and Watson (2012) and Mertens and Ravn
% (2013). 
% =======================================================================
% [IRF, IVout] = VARiriv(VAR,VARopt,IV)
% -----------------------------------------------------------------------
% INPUTS 
%   - VAR: structure, result of VARmodel function
%   - VARopt: options of the VAR (see VARopt from VARmodel)
%   - IV is a (vector) instrumental variable correlated with the shock of
%       interest and uncorrelated with the other shocks
% ----------------------------------------------------------------------- 
% OUTPUT
%   - IRF(t,j): matrix with 't' steps, containing the IRF of 'j' variable 
%       to one identified shock
%   - IVout: first-stage and second stage results from the IV
%       identification
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015; Andrej Sokol, December 2015.
% ambrogio.cesabianchi@gmail.com



%% Check inputs
%===============================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end
if ~exist('IV','var')
    error('You need to provide the instrument (IV)');
end

%% Retrieve and initialize variables 
%=============================================================
nsteps = VARopt.nsteps;
shut   = VARopt.shut;


%% Instrumental variable - First stage
%=========================================================
nvar = VAR.nvar;

% Recover residuals (order matter!)
up = VAR.residuals(:,1);     % is the shock that I want to identify
uq = VAR.residuals(:,2:end); % other shocks

% Make sample of IV comparable with up and uq
[aux, foIV, lo] = CommonSample([up IV(VAR.nlag+1:end)]);
p = aux(:,1);
q = uq(end-length(p)+1:end,:);
Z = aux(:,2:end);

% Run first stage regression and fitted
OLS = OLSmodel(p,Z,1);
p_hat = OLS.yhat;


%% IRF on impact - Second Stage
%=========================================================
IRF(1,1) = 1;
sqsp = zeros(size(q,2),1);
for ii=2:nvar
    second = OLSmodel(q(:,ii-1),p_hat,0);
    IRF(ii,1) = second.beta;
    sqsp(ii-1) = second.beta;
end
if shut~=0
    IRF(shut,1) = 0;
end


%% IRF dynamics
%=========================================================
% Retrieve companion
F = VAR.Fcomp;
if shut~=0
    F(shut,:) = 0;
end

% Make IRF comparable with companion
IRF(end+1:length(F)) = 0;

% Iterate forward
for ii=2:nsteps
    IRF(:,ii) = F*IRF(:,ii-1);
end

% Get IRFiv from big_IRFiv
IRF = IRF(1:nvar,:)';


%% Save
%=========================================================
IVout.first = OLS;
aux = [nan(VAR.nlag,1); nan(foIV,1); p_hat];
IVout.fitted = aux;
IVout.F = OLS.F;
IVout.R2 = OLS.rsqr;

% compute s (Algebra as in Gertler and Karadi, 2015, p. 52)
if ~isempty(varargin)
    sigmaHat = varargin{1};
else
    sigmaHat = VAR.sigma;
end
s21s11 = sqsp; 
S11 = sigmaHat(1,1);
S21 = sigmaHat(2:end,1);
S22 = sigmaHat(2:end,2:end);
Q = s21s11*S11*s21s11'-(S21*s21s11'+s21s11*S21')+S22;
sp = sqrt(S11-(S21-s21s11*S11)'*(Q\(S21-s21s11*S11)));
s = IRF(1,:)'*sp;

IVout.s = s;
IVout.foIV = foIV;



