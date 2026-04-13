%%US VAR
% Set endogenous US
vNamesLongUS = {'1yr Rate';'Bond Yield';'Real GDP';'CPI';'Excess Bond Premium'};
vNamesUS      = {'i_1YR';'BondYield';'lnRGDP';'CPI';'EBP'};
nVarUS   = length(vNamesUS);
vTreatUS = zeros(nVarUS,1);  %0=nothing, 1=log, 2=logdiff, 3=diff
EBPidxUS = find(strcmp(vNamesUS,'EBP'));

% Set IV
vNamesLongIV = {'JarocinskiKaradiMP';'JarocinskiKaradiCBI'};
vNamesIV      = {'MPshockSign';'CBIshockSign'};
vTreatIV      = zeros(length(vNamesIV),1);  %0=nothing, 1=log, 2=logdiff, 3=diff
nVarIV       = length(vNamesIV);

% Load US data
[xlsdata, xlstext] = xlsread(fname,'US');
auxData = Num2NaN(xlsdata);
auxVNames = xlstext(1,2:end);

for ii=1:length(auxVNames)
    
    dataUS.(auxVNames{ii}) = auxData(:,ii);
    
end

date = xlstext(2:end,1);

% Load US IV data
[xlsdata, xlstext] = xlsread(fname,'IV');
auxData = Num2NaN(xlsdata);
auxNames = xlstext(1,2:end);

for ii=1:length(auxNames)
    
    dataIV.(auxNames{ii}) = auxData(:,ii);
    
end

% Observations
fo = find(strcmp(userFO,date));
lo = find(strcmp(userLO,date));
obs = lo-fo+1;

% Create matrices of variables for the VAR and instrument
endoUS = cookVARData(dataUS,obs,nVarUS,fo,lo,1,vTreatUS,vNamesUS);
IV = cookVARData(dataIV,obs,nVarIV,fo,lo,1,vTreatIV,vNamesIV);
VIX = cookVARData(dataUS,obs,1,fo,lo,1,vTreatUS,{'VIX'});

% Estimate VAR and IV IRFs
[varUS, varOptUS] = VARmodel(endoUS,nLag,VARCase);
varOptUS.nsteps = nSteps;
varOptUS.pctg = pctg;
varOptUS.ndraws = nDraws;
varOptUS.method = 'wild';

% External instrument. this pass only needed to work out common sample of
% IV and VAR data
ivUS = VARiriv_B(varUS,IV);
foIV = ivUS.foIV;

% compute sigmaHat using IV sample only
sigmaHat = (1/(varUS.nobs-foIV-varUS.ntotcoeff))*...
    (varUS.residuals(foIV+1:end,:)-repmat(mean(varUS.residuals(foIV+1:end,:)),size(varUS.residuals(foIV+1:end,:),1),1))'*...
    (varUS.residuals(foIV+1:end,:)-repmat(mean(varUS.residuals(foIV+1:end,:)),size(varUS.residuals(foIV+1:end,:),1),1));

% recompute ivUS.s for sigmaHat
ivUS = VARiriv_B(varUS,IV,sigmaHat);

% Sign restrictions
load newSignRestMat

%reshuffle bondyield second
signRestMat = signRestMat([1 5 3 2 4],:); 

% allow more shocks identified with IV
if length(vNamesIV)>1

    signRestMat = signRestMat(:,1:end-length(vNamesIV)+1);
    
end

% variant for unemployment
if strcmp(vNamesUS{3},'U')
    
    signRestMat(3,:) = -signRestMat(3,:);
    
end

if finShockSignOnly
    
    %identify just 1 shock with sign restrictions
    signRestMat(:,2:end) = 0;
    
end

% impose sign restrictions for a longer horizon
signRestMat = repmat(signRestMat,nPeriods,1);

%identification
tic;
[bstUS.draws,bstUS.Fcomp,bstUS.yArtificial,bstUS.srIRF,bstUS.resids,bstUS.A0Mats,...
    bstUS.MTidx,bstUS.nRot, bstUS.R2first,bstUS.relEig] = ...
    VARirbandivSR(varUS,varOptUS,IV,foIV,signRestMat,nPeriods,maxIters);
toc
succDraws = sum(bstUS.nRot>0);
totalRots = sum(bstUS.nRot);
disp(['Bootstraps with at least 1 rotation: ' num2str(succDraws)]);
disp(['Total number of rotations: ' num2str(totalRots)]);
nDraws = length(bstUS.MTidx);

% Compute mean, median and credible intervals
bstUS.inf = prctile(bstUS.srIRF, pctg_inf,4);
bstUS.sup = prctile(bstUS.srIRF, pctg_sup,4);
bstUS.inf2 = prctile(bstUS.srIRF, pctg_inf2,4);
bstUS.sup2 = prctile(bstUS.srIRF, pctg_sup2,4);
bstUS.med = prctile(bstUS.srIRF,50,4);
bstUS.bar = mean(bstUS.srIRF,4);

