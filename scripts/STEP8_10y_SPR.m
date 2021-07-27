%% Long rate differential chart (Figure A1b)
sNames = {'Monetary Policy';'CB Information';'Financial';'Demand';'Supply';'Resid'};
nShocks = 6;
nDraws=1e5; %reduce number as many draws satisfy restrictions
maxIters = 1e2;
% Set endogenous US
vNamesLongUS = {'1yr Rate';'Bond Yield';'Real GDP';'CPI';'Excess Bond Premium';'Long rate'};
vNamesUS      = {'i_1YR';'BondYield';'lnRGDP';'CPI';'EBP';'i_10YR'};
nVarUS   = length(vNamesUS);
vTreatUS = zeros(nVarUS,1);  %0=nothing, 1=log, 2=logdiff, 3=diff

% Set endogenous UK
vNamesLongSOE = {'1yr rate';'CPI';'Exch. Rate ($ per £)';'Corporate Spread';'Real GDP';'Long rate'};
vNamesSOE      = {'i_1YR';'CPI'; 'FX';'CS_i_10YR'; 'ONS_GDP';'i_10YR'};
nVarSOE   = length(vNamesSOE);
vTreatSOE = zeros(nVarSOE,1);  %0=nothing, 1=log, 2=logdiff, 3=diff


%% US VAR
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

% add row and column
signRestMat = [signRestMat zeros(size(signRestMat,1),1)];
signRestMat = [signRestMat;zeros(1,size(signRestMat,2))];
% impose sign restrictions for a longer horizon
signRestMat = repmat(signRestMat,nPeriods,1);

% Identification
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


%% SOE VAR 
[xlsdata, xlstext] = xlsread(fname,'UK');
auxData = Num2NaN(xlsdata);
auxVNames = xlstext(1,2:end);

for ii=1:length(auxVNames)
    
    dataSOE.(auxVNames{ii}) = auxData(:,ii);
    
end

% Create matrices of variables for the VAR
endoSOE = cookVARData(dataSOE,obs,nVarSOE,fo,lo,1,vTreatSOE,vNamesSOE);
starSOE = L(endoUS,-1);

% Define a common sample between all series
[aux, foSOE, ~] = CommonSample([endoSOE starSOE]);

% Define ENDO and EXOG
endoSOE = aux(:,1:nVarSOE);
starSOE = aux(:,nVarSOE+1:end);

% Estimate VAR
[varSOE, varOptSOE] = VARmodel(endoSOE,nLagSOE,VARCase,starSOE,nLagXSOE);

% Construct the full companion.
comp1 = varUS.Fcomp; % north-west block
[rr1, cc1] = size(comp1);
comp4 = varSOE.Fcomp; % south-east block
[rr2, cc2] = size(comp4);
comp2 = zeros(rr1,cc2); % north-east block
exocoeff = varSOE.Ft';
exocoeff = exocoeff(:,VARCase+nVarSOE*nLagSOE+1:end);
[rr3, cc3] = size(exocoeff);
exocoeff = [exocoeff zeros(rr3,cc1-cc3)];
aux = zeros(rr2-rr3,cc1);

if dynamics % south-west block
    
    comp3 = [exocoeff; aux];
    
else
    
    comp3 = [exocoeff; aux];
    comp3 = zeros(size(comp3,1),size(comp3,2));
    
end

% Save the companion and eigenvalues
irfSOE.Fcomp = [comp1 comp2; comp3 comp4];
irfSOE.Eig = max(abs(eig(irfSOE.Fcomp)));

% Initialize IRFs
nTot = sum(bstUS.nRot);
irf = nan(nSteps,nVarSOE,nShocks,nTot);
irfcomp = nan(nVarUS*nLag+nVarSOE*nLagSOE,nShocks,nSteps);
Eig = zeros(nDraws, 1);
bVIXbst = zeros(nDraws,1);

if foSOE>nLag-nLagSOE
    
    startSOE = 1;
    startUS = foSOE-(nLag-nLagSOE)+1;
    
else
    
    startSOE = nLag-nLagSOE-foSOE+1;
    startUS = 1;
    
end

for tt=1:nDraws
    
    % recover US bootstraped draws and generate UK ones
    xArtificial = L(bstUS.yArtificial(:,:,tt),-1);
    xArtificial = xArtificial(foSOE+1:end,:);
    varSOEDraw = wild_bootstrap_VAR(varSOE,startSOE,xArtificial,bstUS.draws(startUS:end,tt));
    
    % Estimate second stage. We do this with the rotations
    % corresponding to each bootstrapped draw from the US
    uUS = squeeze(bstUS.resids(foIV+1:end,:,tt))-...
        repmat(mean(squeeze(bstUS.resids(foIV+1:end,:,tt))),size(squeeze(bstUS.resids(foIV+1:end,:,tt)),1),1);
    
    for mm = 1:bstUS.nRot(tt)
        
        eUS = uUS/(squeeze(bstUS.A0Mats(:,:,sum(bstUS.nRot(1:tt-1))+mm)))';
        uSOE = varSOEDraw.residuals;
        
        if foSOE+nLagSOE>foIV+nLag
            
            aux = CommonSample([uSOE eUS(foSOE+nLagSOE-foIV-nLag+1:end,:)]);
            
        else
            
            aux = CommonSample([uSOE(end-size(eUS,1)+1:end,:) eUS]);
            
        end
        
        uSOE = aux(:,1:nVarSOE);
        fitted = aux(:,nVarSOE+1:end);
        irf(1,:,:,sum(bstUS.nRot(1:tt-1))+mm) = (fitted\uSOE)';
        
        % Construct irfSOE companion on impact
        irfcomp(:,:,1) = ...
            [squeeze(bstUS.srIRF(1,:,:,sum(bstUS.nRot(1:tt-1))+mm)); ...
            zeros(nVarUS*(nLag-1),nShocks); ...
            permute(irf(1,:,:,sum(bstUS.nRot(1:tt-1))+mm),[2 3 1 4]);...
            zeros(nVarSOE*(nLagSOE-1),nShocks)];
        
        % Construct the companion.
        comp1 = bstUS.Fcomp(:,:,tt); % north-west block
        [rr1, cc1] = size(comp1);
        comp4 = varSOEDraw.Fcomp; % south-east block
        [rr2, cc2] = size(comp4);
        comp2 = zeros(rr1,cc2); % north-east block
        exocoeff = varSOEDraw.Ft';
        exocoeff = exocoeff(:,VARCase+nVarSOE*nLagSOE+1:end);
        [rr3, cc3] = size(exocoeff);
        exocoeff = [exocoeff zeros(rr3,cc1-cc3)];
        aux = zeros(rr2-rr3,cc1);
        
        if dynamics
            
            comp3 = [exocoeff; aux];
            
        else
            
            comp3 = [exocoeff; aux];
            comp3 = zeros(size(comp3,1),size(comp3,2));
            
        end
        
        % Save the companion and eigenvalues
        Fcomp = [comp1 comp2; comp3 comp4];
        Eig(tt) = max(abs(eig(Fcomp(:,:))));
        
        if Eig(tt,1)<VARtol
            
            % Compute irfcompSOE for n>1
            for nn=2:nSteps
                
                irfcomp(:,:,nn) = Fcomp(:,:)*irfcomp(:,:,nn-1);
                
            end
            
                irf(:,:,:,sum(bstUS.nRot(1:tt-1))+mm) = permute(irfcomp(nVarUS*nLag+1:nVarUS*nLag+nVarSOE,:,:),...
        [3 1 2]);

            
        end
        
    end
    
end

bstSOE.irf = irf;

% check sign of info shock
for mm = 1:nTot
    
    if bstUS.srIRF(1,3,2,mm)<0 && bstUS.srIRF(1,end-1,2,mm)>0
        
        bstUS.srIRF(:,:,2,mm) = -bstUS.srIRF(:,:,2,mm);
        bstSOE.irf(:,:,2,mm) = -bstSOE.irf(:,:,2,mm);
        
    end

end

% Compute mean, median and credible intervals
longDiff = squeeze(bstSOE.irf(:,strcmp(vNamesSOE,'i_10YR'),:,:) - bstUS.srIRF(:,strcmp(vNamesUS,'i_10YR'),:,:));
bstSOE.longDiffInf = prctile(longDiff,pctg_inf,3);
bstSOE.longDiffSup = prctile(longDiff,pctg_sup,3);
bstSOE.longDiffInf2 = prctile(longDiff,pctg_inf2,3);
bstSOE.longDiffSup2 = prctile(longDiff,pctg_sup2,3);
bstSOE.longDiffMed = prctile(longDiff,50,3);
bstSOE.longDiffBar = nanmean(longDiff,3);


%% Plot
lNames = {'Percent';'Percent';'Percent';'Percent';'Basis Points';'Percent'};
hz = 36;
FigSize(20,14)

for ss = 1
    
    
    H = PlotSwathe(10000*bstSOE.longDiffBar(1:hz,ss), [10000*bstSOE.longDiffSup(1:hz,ss) 10000*bstSOE.longDiffInf(1:hz,ss)], rgb('light red')); hold on;
    H2 = PlotSwathe(10000*bstSOE.longDiffBar(1:hz,ss), [10000*bstSOE.longDiffSup2(1:hz,ss) 10000*bstSOE.longDiffInf2(1:hz,ss)], rgb('dark red'));
    handle(1)=H2.bar; handle(2) = H2.patch; handle(3)=H.patch;
    plot(zeros(hz),'LineWidth',0.5,'LineStyle','-','Color',rgb('black')); hold on
    title(vNamesLongSOE{end})
    
    if ss == 1%3
        ylabel('Basis Points')
    end
    
    xlabel('Months')
    axis tight
    set(gca,'Layer','top')
    optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)
    
end

opt=LegOption;  opt.handle = handle; LegSubplot({'Mean','68% C.I.','90% C. I.'},opt);
optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)
SaveFigure([folderOut 'spread_10y'],highQuality)
close all
