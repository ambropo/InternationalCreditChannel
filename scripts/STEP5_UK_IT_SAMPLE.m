%% UK VAR
[xlsdata, xlstext] = xlsread(fname,'UK');
auxData = Num2NaN(xlsdata);
auxVNames = xlstext(1,2:end);

for ii=1:length(auxVNames)
    
    dataSOE.(auxVNames{ii}) = auxData(:,ii);
    
end

% Create matrices of variables for the VAR
endoSOE = cookVARData(dataSOE,obs,nVarSOE,fo,lo,1,vTreatSOE,vNamesSOE);
starSOE = L(endoUS,-1);

%IT sample, starts in 1993
endoSOE(1:end-24*12,:) = NaN; 

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
irfcomp = nan(nVarUS*nLag+nVarSOE*nLagSOE,nShocks,nSteps,nTot);
Fcomp = zeros([size(irfSOE.Fcomp) nDraws]);
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
        irfcomp(:,:,1,sum(bstUS.nRot(1:tt-1))+mm) = ...
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
        Fcomp(:,:,tt) = [comp1 comp2; comp3 comp4];
        Eig(tt) = max(abs(eig(Fcomp(:,:,tt))));
        
        if Eig(tt,1)<VARtol
            
            % Compute irfcompSOE for n>1
            for nn=2:nSteps
                
                irfcomp(:,:,nn,sum(bstUS.nRot(1:tt-1))+mm) = Fcomp(:,:,tt)*irfcomp(:,:,nn-1,sum(bstUS.nRot(1:tt-1))+mm);
                
            end
            
        end
        
    end
    
end

for mm = 1:nTot
    
    irf(:,:,:,mm) = permute(squeeze(irfcomp(nVarUS*nLag+1:nVarUS*nLag+nVarSOE,:,:,mm)),...
        [3 1 2]);
    
end

bstSOE.irf = irf;

% check sign of info shock
for mm = 1:nTot
    
    if bstUS.srIRF(1,3,2,mm)<0 && bstUS.srIRF(1,end,2,mm)>0
        
        bstUS.srIRF(:,:,2,mm) = -bstUS.srIRF(:,:,2,mm);
        bstSOE.irf(:,:,2,mm) = -bstSOE.irf(:,:,2,mm);
        
    end

end

% Compute mean, median and credible intervals
bstSOE.inf = prctile(bstSOE.irf,pctg_inf,4);
bstSOE.sup = prctile(bstSOE.irf,pctg_sup,4);
bstSOE.inf2 = prctile(bstSOE.irf,pctg_inf2,4);
bstSOE.sup2 = prctile(bstSOE.irf,pctg_sup2,4);
bstSOE.med = prctile(bstSOE.irf,50,4);
bstSOE.bar = nanmean(bstSOE.irf,4);


%% Plots
lNames = {'Percent';'Percent';'Percent';'Percent';'Basis Points'};
hz = 36;

% UK shocks
FigSize(28,20)

for ss = 1:3
    
    index = 1;
    
    for hh=1:nVarSOE
        
        subplot(2,3,index)
               
        H = PlotSwathe(100*bstSOE.bar(1:hz,hh,ss), [100*bstSOE.sup(1:hz,hh,ss) 100*bstSOE.inf(1:hz,hh,ss)], rgb('light red')); hold on;
        H2 = PlotSwathe(100*bstSOE.bar(1:hz,hh,ss), [100*bstSOE.sup2(1:hz,hh,ss) 100*bstSOE.inf2(1:hz,hh,ss)], rgb('dark red'));
        handle(1)=H2.bar; handle(2) = H2.patch; handle(3)=H.patch;
        plot(zeros(nSteps),'LineWidth',0.5,'LineStyle','-','Color',rgb('black')); hold on
        title(vNamesLongSOE(hh))
        
        if hh == 4
            
            ylabel('Basis Points')
            
        else
            
            ylabel('Percent')
            
        end
        
        xlabel('Months')
        axis tight
        set(gca,'Layer','top')
        optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)
        index = index+1;
        
    end
    
    opt=LegOption; opt.handle = handle; LegSubplot({'Mean','68% C.I.','90% C. I.'},opt);
    optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)
    SaveFigure([folderOut 'UK ' sNames{ss} ' IT sample'],highQuality)
    clf('reset')
    
end

close all
