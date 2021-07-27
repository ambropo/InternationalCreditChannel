function VARout = wild_bootstrap_VAR(VARin,startSOE,xArtificial,draws)
% spits out a VAR estimated based on bootstrapped residuals given by
% "draws", which is a vector of +1 and -1, following wild bootstrap
% procedure.

Ft      = VARin.Ft;  % rows are coefficients, columns are equations
nVar    = VARin.nvar;
nVarExo = VARin.nvar_ex;
nLagXSOE = VARin.nlag_ex;
nLagSOE    = VARin.nlag;
const   = VARin.const;
resid   = VARin.residuals;
endo    = VARin.ENDO;

%% STEP 1: bootstrapped residuals (consistent with US draw)
u = resid(startSOE:end,:).*(draws*ones(1,nVar));
nObs = size(u,1);
%% STEP 2.1: initial values for the artificial data
% Intialize the first nlag observations with real data
LAG=[];
LAGX = [];
for jj = startSOE:startSOE-1+nLagSOE
    yArtificial(jj,:) = endo(jj,:);
    LAG = [yArtificial(jj,:) LAG];
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
   
    LAGplus = [1 T(1) T(1).^2 LAG];
end

if nVarExo~=0
    for jj = startSOE+nLagSOE-nLagXSOE:startSOE+nLagSOE
        LAGX = [xArtificial(jj,:) LAGX];
    end
    LAGplus = [LAGplus LAGX]; 
end

%% STEP 2.2: generate artificial series
% From observation nlag+1 to nobs, compute the artificial data
for jj = nLagSOE+startSOE:nObs+nLagSOE
    for mm = 1:nVar
        % Compute the value for time=jj
        yArtificial(jj,mm) = LAGplus * Ft(1:end,mm) + u(jj-nLagSOE,mm);
    end
    % now update the LAG matrix
    if jj<nObs+nLagSOE
        LAG = [yArtificial(jj,:) LAG(1,1:(nLagSOE-1)*nVar)];
               
        if const==0
            LAGplus = LAG;
        elseif const==1
            LAGplus = [1 LAG];
        elseif const==2
            LAGplus = [1 T(jj-nLagSOE+1) LAG];
        elseif const==3
            LAGplus = [1 T(jj-nLagSOE+1) T(jj-nLagSOE+1).^2 LAG];
        end
        if nVarExo~=0
            if nLagXSOE>0
                LAGX = [xArtificial(jj,:) LAGX(1,1:nLagXSOE*nVarExo)];
            else
                LAGX = xArtificial(jj,:);
            end
            LAGplus = [LAGplus LAGX];
        end
    end
end

%% STEP 3: estimate VAR on artificial data.
if nVarExo~=0
    [VARout, ~] = VARmodel(yArtificial,nLagSOE,const,xArtificial(startSOE:end,:),nLagXSOE);
else
    [VARout, ~] = VARmodel(yArtificial,nLagSOE,const);
end

end