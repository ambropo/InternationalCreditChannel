function [RR,Fcomp,yArtificial,srIRF,resids,...
    A0Mats,MTidx,nRot,R2first,relEig,Ffirst] = VARirbandivSR(VAR,VARopt,IV,foIV,signRestMat,...
    nPeriods,maxIters,varargin)
% =======================================================================
% Calculates confidence intervals for impulse response functions computed
% with VARiriv
% =======================================================================
% [INF,SUP,MED,BAR] = VARirbandiv(VAR,VARopt,IV)
% -----------------------------------------------------------------------
% INPUTS 
%   - VAR   : VAR results obtained with VARmodel (structure)
%	- VARopt: options of the IRFs (see VARoption)
% -----------------------------------------------------------------------
% OUTPUT
%   - INF(t,j): lower confidence band (t steps, j variable)
%   - SUP(t,j): upper confidence band (t steps, j variable)
%   - MED(t,j): median response (t steps, j variable)
%   - BAR(t,j): mean response (t steps, j variable)
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com


%% Check inputs
%===============================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end
if ~exist('IV','var')
    error('You need to provide an instrumental variable (IV)');
end


%% Retrieve and initialize variables 
%=============================================================
% nsteps = VARopt.nsteps;
nDraws = VARopt.ndraws;
method = VARopt.method;
nI = size(IV,2);

Ft      = VAR.Ft;  % rows are coefficients, columns are equations
nVar    = VAR.nvar;
nVar_ex = VAR.nvar_ex;
nLag    = VAR.nlag;
const   = VAR.const;
nObs    = VAR.nobs;
resid   = VAR.residuals;
ENDO    = VAR.ENDO;
EXOG    = VAR.EXOG;
nTotCoeff = VAR.ntotcoeff;

srIRF= cell(nDraws,1);%double.empty(nsteps,nVar,nVar,0);
MTidx = zeros(nDraws,1);
nRot = MTidx;
resids = cell(nDraws,1);%zeros([size(resid) nDraws]);
R2first = [];%NaN(nDraws,nI);
relEig = R2first;
Ffirst = R2first;

%% Create the matrices for the loop
%==================================
% y_artificial = zeros(nObs+nLag,nVar);
yArtificial = cell(nDraws,1);%NaN([size(y_artificial) nDraws]);
RR = NaN(nObs,nDraws);
Fcomp = cell(nDraws,1);%NaN(nVar*nLag,nVar*nLag,nDraws);

%% Loop over the number of draws, generate data, estimate the var and then 
%% calculate impulse response functions
%==========================================================================
% ww = 1; % index for printing on screen
disp('Starting draws...')
parfor tt=1:nDraws

    y_artificial = zeros(nObs+nLag,nVar);
    
%     %Display number of loops
%     if tt==100*ww
%         disp(['Loop ' num2str(tt) ' / ' num2str(nDraws) ' draws'])
%         ww=ww+1;
%     end

%% STEP 1: choose the method (Monte Carlo or Bootstrap) and generate the 
%% residuals
    if strcmp(method,'bs')
        % Use the residuals to bootstrap: generate a random number bounded 
        % between 0 and # of residuals, then use the ceil function to select 
        % that row of the residuals (this is equivalent to sampling with replacement)
        error(['The method ' method ' is not available for VARirbandiv. Use wild boostrap.'])
    elseif strcmp(method,'wild')
        % Wild bootstrap based on simple distribution (~Rademacher)
        rr = 1-2*(rand(nObs,1)>0.5);
        u = resid.*(rr*ones(1,nVar));
        Z = [IV(1:nLag,:); IV(nLag+1:end,:).*rr];
    else
        error(['The method ' method ' is not available'])
    end

%% STEP 2: generate the artifcial data

    %% STEP 2.1: initial values for the artificial data
    % Intialize the first nlag observations with real data
    LAG=[];
    for jj = 1:nLag
        y_artificial(jj,:) = ENDO(jj,:);
        LAG = [y_artificial(jj,:) LAG]; 
    end
    % Initialize the artificial series and the LAGplus vector
    T = [1:nObs]';
    if const==0
        LAGplus = LAG;
    elseif const==1
        LAGplus = [1 LAG];
    elseif const==2
        LAGplus = [1 T(1) LAG]; 
    elseif const==3
        T = [1:nObs]';
        LAGplus = [1 T(1) T(1).^2 LAG];
    end
    if nVar_ex~=0
        LAGplus = [LAGplus EXOG(jj+1,:)]; %not sure the timing is right here, though it doesn;t matter
    end
    
    %% STEP 2.2: generate artificial series
    % From observation nlag+1 to nobs, compute the artificial data
    for jj = nLag+1:nObs+nLag
        for mm = 1:nVar
            % Compute the value for time=jj
            y_artificial(jj,mm) = LAGplus * Ft(1:end,mm) + u(jj-nLag,mm);
        end
        % now update the LAG matrix
        if jj<nObs+nLag
            LAG = [y_artificial(jj,:) LAG(1,1:(nLag-1)*nVar)];
            if const==0
                LAGplus = LAG;
            elseif const==1
                LAGplus = [1 LAG];
            elseif const==2
                LAGplus = [1 T(jj-nLag+1) LAG];
            elseif const==3
                LAGplus = [1 T(jj-nLag+1) T(jj-nLag+1).^2 LAG];
            end
            if nVar_ex~=0
                LAGplus = [LAGplus EXOG(jj+1,:)];
            end
        end
    end

%% STEP 3: estimate VAR on artificial data. 
    if nVar_ex~=0
        [VAR_draw, ~] = VARmodel(y_artificial,nLag,const,EXOG);
    else
        [VAR_draw, ~] = VARmodel(y_artificial,nLag,const);
    end
    
%% STEP 4: calculate "ndraws" impulse responses and store them

%     IV_draw = VARiriv_B(VAR_draw,Z);  % uses options from VARopt, but companion etc. from VAR_draw
    % should initialise all these variables....
    if VAR_draw.maxEig<.999
        
       RR(:,tt) = rr;
       sigmaHat = (1/(nObs-foIV-nTotCoeff))*(VAR_draw.residuals(foIV+1:end,:)-repmat(mean(VAR_draw.residuals(foIV+1:end,:)),size(VAR_draw.residuals(foIV+1:end,:),1),1))'*...
       (VAR_draw.residuals(foIV+1:end,:)-repmat(mean(VAR_draw.residuals(foIV+1:end,:)),size(VAR_draw.residuals(foIV+1:end,:),1),1));
%         sigmaHat = VAR_draw.sigma;
        IV_draw = VARiriv_B(VAR_draw,Z,sigmaHat);  % uses options from VARopt, but companion etc. from VAR_draw
       [auxIRF,~,~,~,~,MTidx(tt)] = VARirSR(VAR_draw,VARopt,IV_draw.s,sigmaHat,signRestMat,nPeriods,maxIters);
%        R2first(tt,:) = IV_draw.R2first;
%        Ffirst(tt,:) = IV_draw.Ffirst;
%        relEig(tt,:) = sort(eig(IV_draw.Lambda),'descend');
       
    end

        if MTidx(tt) > 0
            
            Fcomp{tt} = VAR_draw.Fcomp;
            resids{tt} = VAR_draw.residuals;
            yArtificial{tt} = y_artificial;
            srIRF{tt}=auxIRF;%cat(4,srIRF,permute(auxIRF,[2 3 4 1]));
            nRot(tt) = size(auxIRF,1);
            
            
        end
        
end

disp('-- Done!');
srIRF = permute(cat(1,srIRF{MTidx>0}),[2 3 4 1]);
Fcomp = cat(3,Fcomp{MTidx>0});
resids = cat(3,resids{MTidx>0});
yArtificial = cat(3,yArtificial{MTidx>0});


% extract impact matrices
A0Mats = squeeze(srIRF(1,:,:,:));
RR = RR(:,MTidx>0);
% Fcomp = Fcomp(:,:,MTidx>0);
% resids = resids(:,:,MTidx>0);
% yArtificial = yArtificial(:,:,MTidx>0);
MTidx = MTidx(MTidx>0);
nRot = nRot(nRot>0);

end
