function [IVout] = VARiriv_B(VAR,IV,varargin)
% =======================================================================
% Compute IRFs for a VAR model estimated with VARmodel using the external 
% instrument approach of Stock and Watson (2012) and Mertens and Ravn
% (2013). 
% =======================================================================
% [IVout] = VARiriv(VAR,VARopt,IV)
% -----------------------------------------------------------------------
% INPUTS 
%   - VAR: structure, result of VARmodel function
%   - VARopt: options of the VAR (see VARopt from VARmodel)
%   - resIdx: picks residuals for first stage regression
%   - IV is a (vector) instrumental variable correlated with the shocks of
%       interest and uncorrelated with the other shocks
% ----------------------------------------------------------------------- 
% OUTPUT
%   - IVout: first-stage and second stage results from the IV
%       identification
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015; Andrej Sokol, October 2020.
% ambrogio.cesabianchi@gmail.com



%% Check inputs
%===============================================
if ~exist('IV','var')
    error('You need to provide the instrument (IV)');
end


%% Instrumental variable - First stage
%=========================================================
nShocks = size(IV,2);

% Recover residuals (order matter!)
up = VAR.residuals(:,1:nShocks);     % is the shock that I want to identify
uq = VAR.residuals(:,nShocks+1:end); % other shocks

% Make sample of IV comparable with up and uq
[aux, foIV] = CommonSample([up IV(VAR.nlag+1:end,:)]);
u1 = aux(:,1:nShocks);
T = size(u1,1);
u2 = uq(end-T+1:end,:);
u1 = u1 - mean(u1); %this is needed because IV sample shorter than VAR sample. IV assumed to be mean 0 here...
u2 = u2 - mean(u2);
m = aux(:,nShocks+1:end);

% Run first & second stage regressions
bOLS = [ones(T,1) m]\u1;%OLSmodel(u1(:,1),m,1);
u1Hat = [ones(T,1) m]*bOLS;%OLS.yhat;
% bOLS = m\u1;%OLSmodel(u1(:,1),m,1);
% u1Hat = m*bOLS;%OLS.yhat;

% bOLS = bOLS(2:end,:);

invSmu1Smu2 = u1Hat\u2;

% compute b(eta) (s in our notation) (Algebra as in Mertens and Ravn, 2013, pp. 1234-5)
if ~isempty(varargin)

    sigmaHat = varargin{1};

else
    
    sigmaHat = VAR.sigma;

end

% elements of varcov matrix (re-index to also allow 1 instr for checks)
S11 = sigmaHat(1:nShocks,1:nShocks);
S21 = sigmaHat(nShocks+1:end,1:nShocks);
S22 = sigmaHat(nShocks+1:end,nShocks+1:end);

% auxiliary terms
b21invb11 = invSmu1Smu2'; %check transposition
b21invb11p = invSmu1Smu2;
Z = (b21invb11*S11)*b21invb11p - (S21*b21invb11p+b21invb11*S21') + S22;
b12b12p = (S21-b21invb11*S11)'*inv2(Z)*(S21-b21invb11*S11);
b22b22p = S22 + (b21invb11*(b12b12p-S11))*b21invb11p;
b12invb22 = (b12b12p*b21invb11p + (S21 - b21invb11*S11)')/b22b22p; %check last term, typo in paper?
b11b11p = S11-b12b12p;

%A1-A3
b11invS1 = inv2(eye(nShocks) - b12invb22*b21invb11);
b21invS1 = b21invb11*inv2(eye(nShocks) - b12invb22*b21invb11);
S1S1p = ((eye(nShocks) - b12invb22*b21invb11)*b11b11p)*(eye(nShocks) - b12invb22*b21invb11)';

try 
    
    S1 = chol(S1S1p,'lower');
    
catch
    
    [V,D] = eig(S1S1p);
    D = diag(D);
    D(D<0) = 1e-10;
    S1 = chol(V*(V'.*D),'lower');

end

b11 = b11invS1*S1;
b21 = b21invS1*S1;
b = [b11;b21];

IVout.s = b; 
IVout.foIV = foIV;

% Instrument reliability statistics (see Mertens and Ravn, 2013, Appendix and Table on p. 1229)
u1Res = u1 - u1Hat;
nonCensored = sum(abs(m),2)>0;
d = sum(nonCensored)/length(nonCensored);
Smm = (m'*m)/T;
Smu1 = (m'*u1)/T;
IVout.Lambda = inv2(Smm)*Smu1*(inv2(b11b11p)*Smu1')/d;
IVout.R2first = diag(1-(u1Res(nonCensored,:)'*u1Res(nonCensored,:))./(u1(nonCensored,:)'*u1(nonCensored,:)));
IVout.Ffirst = diag(((u1(nonCensored,:)'*u1(nonCensored,:)-u1Res(nonCensored,:)'*u1Res(nonCensored,:))/nShocks)./((u1Res(nonCensored,:)'*u1Res(nonCensored,:)/(sum(nonCensored)-nShocks-1))));
% IVout.Ffirst = diag(((u1'*u1-u1Res'*u1Res)/nShocks)./((u1Res'*u1Res/(length(u1)-nShocks-1))));



