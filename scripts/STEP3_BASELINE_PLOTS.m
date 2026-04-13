%% Baseline Plots (Figures 4-9 and A1a)
lNames = {'Percent';'Percent';'Percent';'Percent';'Basis Points'};
hz = 36;

% US shocks
FigSize(28,20)

for ss = 1:3
    
    index = 1;
    
    for hh=1:nVarUS
        
        subplot(2,3,index)        
        H = PlotSwathe(100*squeeze(bstUS.bar(1:hz,hh,ss)), [100*bstUS.sup(1:hz,hh,ss) 100*bstUS.inf(1:hz,hh,ss)], rgb('light blue')); hold on;
        H2 = PlotSwathe(100*squeeze(bstUS.bar(1:hz,hh,ss)), [100*bstUS.sup2(1:hz,hh,ss) 100*bstUS.inf2(1:hz,hh,ss)], rgb('dark blue'));
        handle(1)=H2.bar; handle(2) = H2.patch; handle(3)=H.patch;
        plot(zeros(nSteps),'LineWidth',0.5,'LineStyle','-','Color',rgb('black')); hold on
        title(vNamesLongUS(hh))
        ylabel(lNames{hh})
        xlabel('Months')
        axis tight
        set(gca,'Layer','top')
        optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)
        index = index+1;
        
    end
    
    opt=LegOption; opt.handle = handle; LegSubplot({'Mean','68% C.I.','90% C. I.'},opt);
    optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)
    SaveFigure([folderOut 'US ' sNames{ss}],highQuality)
    clf('reset')
    
end

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
    SaveFigure([folderOut 'UK ' sNames{ss}],highQuality)
    clf('reset')
    
end

close all

% spreads panel
FigSize(20,14)

for ss = 1
    
    H = PlotSwathe(10000*bstSOE.shortDiffBar(1:hz,ss), [10000*bstSOE.shortDiffSup(1:hz,ss) 10000*bstSOE.shortDiffInf(1:hz,ss)], rgb('light red')); hold on;
    H2 = PlotSwathe(10000*bstSOE.shortDiffBar(1:hz,ss), [10000*bstSOE.shortDiffSup2(1:hz,ss) 10000*bstSOE.shortDiffInf2(1:hz,ss)], rgb('dark red'));    
    handle(1)=H2.bar; handle(2) = H2.patch; handle(3)=H.patch;
    plot(zeros(hz),'LineWidth',0.5,'LineStyle','-','Color',rgb('black')); hold on
    title(vNamesLongSOE{1})
    
    if ss == 1
        ylabel('Basis Points')
    end
    
    xlabel('Months')
    axis tight
    set(gca,'Layer','top')
    optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)
    
end

opt=LegOption;  opt.handle = handle; LegSubplot({'Mean','68% C.I.','90% C. I.'},opt);
optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)

SaveFigure([folderOut 'spread_1y'],highQuality)
close all

