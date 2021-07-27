%% UK VAR (BASELINE AND COUNTERFACTUAL)

% Set endogenous UK
vNamesLongSOE = {'1yr rate';'CPI';'Exch. Rate ($ per £)';'Corporate Spread';'Real GDP'};
vNamesSOE      = {'i_1YR';'CPI'; 'FX';'CS_i_10YR'; 'ONS_GDP'};
nVarSOE   = length(vNamesSOE);
vTreatSOE = zeros(nVarSOE,1);  %0=nothing, 1=log, 2=logdiff, 3=diff
EBPidxSOE = find(strcmp(vNamesSOE,'CS_i_10YR'));
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
    bstSOE.irf = nan(nSteps,nVarSOE,nShocks,nTot);
    bstSOE.irfCF = bstSOE.irf;
    bstUS.irfCF = nan(nSteps,nVarUS,nShocks,nTot);
    bstSOE.irfcomp = nan(nVarUS*nLag+nVarSOE*nLagSOE,nShocks,nSteps,nTot);
    bstSOE.irfcompCF = bstSOE.irfcomp;
    bstSOE.Fcomp = zeros([size(irfSOE.Fcomp) nDraws]);
    bstSOE.FcompCF = bstSOE.Fcomp;
    bstSOE.Eig = zeros(nDraws, 1);
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
            bstSOE.irf(1,:,:,sum(bstUS.nRot(1:tt-1))+mm) = (fitted\uSOE)';
            
            % Construct irfSOE companion on impact
            bstSOE.irfcomp(:,:,1,sum(bstUS.nRot(1:tt-1))+mm) = ...
                [squeeze(bstUS.srIRF(1,:,:,sum(bstUS.nRot(1:tt-1))+mm)); ...
                zeros(nVarUS*(nLag-1),nShocks); ...
                permute(bstSOE.irf(1,:,:,sum(bstUS.nRot(1:tt-1))+mm),[2 3 1 4]);...
                zeros(nVarSOE*(nLagSOE-1),nShocks)];
            bstSOE.irfcompCF(:,:,1,sum(bstUS.nRot(1:tt-1))+mm) = bstSOE.irfcomp(:,:,1,sum(bstUS.nRot(1:tt-1))+mm);
            % CF response in both countries on impact
            bstSOE.irfcompCF([EBPidxUS nVarUS*nLag+EBPidxSOE],:,1,sum(bstUS.nRot(1:tt-1))+mm) = 0;
            
            % Construct the companion. 
            comp1 = bstUS.Fcomp(:,:,tt); % north-west block
            comp1CF = comp1;
            comp1CF(EBPidxUS,:) = 0;
            [rr1, cc1] = size(comp1);
            comp4 = varSOEDraw.Fcomp; % south-east block
            comp4CF = comp4;
            comp4CF(EBPidxSOE,:) = 0;
            [rr2, cc2] = size(comp4);
            comp2 = zeros(rr1,cc2); % north-east block
            exocoeff = varSOEDraw.Ft';
            exocoeff = exocoeff(:,VARCase+nVarSOE*nLagSOE+1:end);
            [rr3, cc3] = size(exocoeff);
            exocoeff = [exocoeff zeros(rr3,cc1-cc3)];
            aux = zeros(rr2-rr3,cc1);
            
            if dynamics
                
                comp3 = [exocoeff; aux];
                comp3CF = comp3;
                comp3CF(EBPidxSOE,:) = 0;

            
            else
                
                comp3 = [exocoeff; aux];
                comp3 = zeros(size(comp3,1),size(comp3,2)); 
                
            end
            
            % Save the companion and eigenvalues
            bstSOE.Fcomp(:,:,tt) = [comp1 comp2; comp3 comp4];
            bstSOE.FcompCF(:,:,tt) = [comp1CF comp2; comp3CF comp4CF];
            bstSOE.Eig(tt) = max(abs(eig(bstSOE.Fcomp(:,:,tt))));
            
            if bstSOE.Eig(tt,1)<VARtol
                
                % Compute irfcompSOE for n>1
                for nn=2:nSteps
                    
                    bstSOE.irfcomp(:,:,nn,sum(bstUS.nRot(1:tt-1))+mm) = bstSOE.Fcomp(:,:,tt)*bstSOE.irfcomp(:,:,nn-1,sum(bstUS.nRot(1:tt-1))+mm);
                    % counterfactual
                    bstSOE.irfcompCF(:,:,nn,sum(bstUS.nRot(1:tt-1))+mm) = bstSOE.FcompCF(:,:,tt)*bstSOE.irfcompCF(:,:,nn-1,sum(bstUS.nRot(1:tt-1))+mm);

                end
                
            end
            
        end

    end
    
        
    for mm = 1:nTot

        bstSOE.irf(:,:,:,mm) = permute(squeeze(bstSOE.irfcomp(nVarUS*nLag+1:nVarUS*nLag+nVarSOE,:,:,mm)),...
            [3 1 2]);
        bstSOE.irfCF(:,:,:,mm) = permute(squeeze(bstSOE.irfcompCF(nVarUS*nLag+1:nVarUS*nLag+nVarSOE,:,:,mm)),...
                [3 1 2]);
            bstUS.irfCF(:,:,:,mm) = permute(squeeze(bstSOE.irfcompCF(1:nVarUS,:,:,mm)),...
                [3 1 2]);

    end
    
% check sign of info shock
for mm = 1:nTot
    
    if bstUS.srIRF(1,3,2,mm)<0 && bstUS.srIRF(1,end,2,mm)>0
        
        bstUS.srIRF(:,:,2,mm) = -bstUS.srIRF(:,:,2,mm);
        bstSOE.irf(:,:,2,mm) = -bstSOE.irf(:,:,2,mm);
        
    end
    
    if bstUS.irfCF(1,3,2,mm)<0 && bstUS.irfCF(1,end,2,mm)>0
        
        bstUS.irfCF(:,:,2,mm) = -bstUS.irfCF(:,:,2,mm);
        bstSOE.irfCF(:,:,2,mm) = -bstSOE.irfCF(:,:,2,mm);
        
    end
    
end

% Compute mean, median and credible intervals
bstSOE.inf = prctile(bstSOE.irf,pctg_inf,4);
bstSOE.sup = prctile(bstSOE.irf,pctg_sup,4);
bstSOE.inf2 = prctile(bstSOE.irf,pctg_inf2,4);
bstSOE.sup2 = prctile(bstSOE.irf,pctg_sup2,4);
bstSOE.med = prctile(bstSOE.irf,50,4);
bstSOE.bar = nanmean(bstSOE.irf,4);
bstSOE.medCF = prctile(bstSOE.irfCF,50,4);
bstSOE.barCF = nanmean(bstSOE.irfCF,4);
bstUS.medCF = prctile(bstUS.irfCF,50,4);
bstUS.barCF = nanmean(bstUS.irfCF,4);
shortDiff = squeeze(bstSOE.irf(:,strcmp(vNamesSOE,'i_1YR'),:,:) - bstUS.srIRF(:,strcmp(vNamesUS,'i_1YR'),:,:));
bstSOE.shortDiffInf = prctile(shortDiff,pctg_inf,3);
bstSOE.shortDiffSup = prctile(shortDiff,pctg_sup,3);
bstSOE.shortDiffInf2 = prctile(shortDiff,pctg_inf2,3);
bstSOE.shortDiffSup2 = prctile(shortDiff,pctg_sup2,3);
bstSOE.shortDiffMed = prctile(shortDiff,50,3);
bstSOE.shortDiffBar = nanmean(shortDiff,3);
