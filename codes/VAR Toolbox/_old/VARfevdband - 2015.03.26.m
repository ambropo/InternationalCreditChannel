function [INF,SUP,MED] = VARfevdband(VAR,VARopt)
% =======================================================================
% Calculate confidence intervals for forecast error variance decomposition
% computed with VARfevd
% =======================================================================
% [INF,SUP,MED] = VARfevdband(VAR,FEVD_opt,ndraws,pctg,method)
% -----------------------------------------------------------------------
% INPUTS 
%   - VAR: VAR results obtained with VARmodel (structure)
%	- VARopt: options of the FEVDs (see VARoption)
% -----------------------------------------------------------------------
% OUTPUT
%   - INF(t,j,k): lower confidence band (t steps, j variable, k shock)
%   - SUP(t,j,k): upper confidence band (t steps, j variable, k shock)
%   - MED(t,j,k): median response (t steps, j variable, k shock)
% =======================================================================
% Ambrogio Cesa Bianchi, march 2015
% ambrogio.cesabianchi@gmail.com

% I thank Fabien Tripier for finding a bug in STEP 2.2 for the case of
% constant and trend (det==2).


%% Check inputs
%===============================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end


%% Retrieve and initialize variables 
%=============================================================
nsteps = VARopt.nsteps;
impact = VARopt.impact;
ident  = VARopt.ident;
ndraws = VARopt.ndraws;
pctg   = VARopt.pctg;
method = VARopt.method;

Ft      = VAR.Ft; % rows are coefficients, columns are equations
nvars   = VAR.nvar;
nvar_ex = VAR.nvar_ex;
nlag    = VAR.nlag;
det     = VAR.det;
nobs    = VAR.nobs;
Y       = VAR.Y;
resid   = VAR.residuals;
if nvar_ex~=0
    exog = VAR.X_EX;
end

INF = zeros(nsteps,nvars,nvars);
SUP = zeros(nsteps,nvars,nvars);
MED = zeros(nsteps,nvars,nvars);

%% Create the matrices for the loop
%==================================
y_artificial = zeros(nobs,nvars);

%% Loop over the number of draws, generate data, estimate the var and then 
%% calculate impulse response functions
%==========================================================================

for tt=1:ndraws
    disp(['Loop ' num2str(tt) ' / ' num2str(ndraws) ' draws'])

%% STEP 1: choose the method (Monte Carlo or Bootstrap) and generate the
%% residuals
    if strcmp(method,'bs')
        % Use the residuals to bootstrap: generate a random number bounded 
        % between 0 and # of residuals, then use the ceil function to select 
        % that row of the residuals (this is equivalent to sampling with replacement)
        u = resid(ceil(size(resid,1)*rand(nobs,1)),:);
    else
        error(['The method ' method ' is not available'])
    end

%% STEP 2: generate the artifcial data

    %% STEP 2.1: generate initial values for the artifcial data
    % Intialize the first nlag observations with real data + plus artificial
    % res. Nontheless, in the estimation of the var on the simulated data, 
    % I through away the first nobs observations so it should not matter.
    LAG=[];
    for jj = 1:nlag
        y_artificial(jj,:) = Y(jj,:) + u(jj,:);
        LAG = [y_artificial(jj,:) LAG]; 
        % Initialize the artificial series and the LAGplus vector
        if det==0
            LAGplus = LAG;
        elseif det==1
            LAGplus = [1 LAG];
        elseif det==2
            T = [1:nobs]';
            LAGplus = [1 T(jj) LAG];
        elseif det==3
            T = [1:nobs]';
            LAGplus = [1 T(jj) T(jj).^2 LAG];
        end
        if nvar_ex~=0
            LAGplus = [LAGplus exog(jj,:)];
        end
    end
    
    %% STEP 2.2: generate artificial series 
    % From observation nlag+1 to nobs, compute the artificial data
    for jj = nlag+1:nobs
        for mm = 1:nvars
            % Compute the value for time=jj
            y_artificial(jj,mm) = LAGplus * Ft(1:end,mm) + u(jj,mm);
        end
        % now update the LAG matrix
        LAG = [y_artificial(jj,:) LAG(1,1:(nlag-1)*nvars)];
        if det==0
            LAGplus = LAG;
        elseif det==1
            LAGplus = [1 LAG];
        elseif det==2
            LAGplus = [1 T(jj) LAG];
        elseif det==3
            LAGplus = [1 T(jj) T(jj).^2 LAG];
        end
        if nvar_ex~=0
            LAGplus = [LAGplus exog(jj,:)];
        end
    end

%% STEP 3: estimate VAR on artificial data
    if nvar_ex~=0
        VAR_draw = VARmodel(y_artificial(1:end,:),nlag,det,exog);
    else
        VAR_draw = VARmodel(y_artificial(1:end,:),nlag,det);
    end
    
    beta_draw  = VAR_draw.Ft;
    sigma_draw = VAR_draw.sigma;

%% STEP 4: calculate "ndraws" fevd and store them

    [fevd_draw, ~] = VARfevd(VAR_draw,VARopt);
    
    % if you don't have three dimensional arrays this will break.
    FEVD(:,:,:,tt) = fevd_draw;
    
end

%% Compute the error bands
%=========================
if strcmp(method,'bs') % When using boostratp, use percentile (upper and lower bounds) bands type
    pctg_inf = (100-pctg)/2; 
    pctg_sup = 100 - (100-pctg)/2;
    INF(:,:,:) = prctile(FEVD(:,:,:,:),pctg_inf,4);
    SUP(:,:,:) = prctile(FEVD(:,:,:,:),pctg_sup,4);
    MED(:,:,:) = prctile(FEVD(:,:,:,:),50,4);
else
    error(['The method ' method ' is not available'])
end

