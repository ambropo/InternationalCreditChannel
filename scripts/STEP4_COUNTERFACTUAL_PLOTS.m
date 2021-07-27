%% Counterfactual Plots (Figures 10 and B4-B9)
lNames = {'Percent';'Percent';'Percent';'Percent';'Basis Points'};
hz = 36;

% Paper plot
FigSize(28,20)
usGDP = strcmp(vNamesUS,'lnRGDP');
ukGDP = strcmp(vNamesSOE,'ONS_GDP');


index = 1;

for cc = 1:2

    for ss = [3 1 2]
                
        subplot(2,3,index)
        
        if cc == 1
            
            H = PlotSwathe(100*squeeze(bstUS.bar(1:hz,usGDP,ss)), [100*bstUS.sup(1:hz,usGDP,ss) 100*bstUS.inf(1:hz,usGDP,ss)], rgb('light blue')); hold on;
            H2 = PlotSwathe(100*squeeze(bstUS.bar(1:hz,usGDP,ss)), [100*bstUS.sup2(1:hz,usGDP,ss) 100*bstUS.inf2(1:hz,usGDP,ss)], rgb('dark blue'));
            
        else
            H = PlotSwathe(100*bstSOE.bar(1:hz,ukGDP,ss), [100*bstSOE.sup(1:hz,ukGDP,ss) 100*bstSOE.inf(1:hz,ukGDP,ss)], rgb('light red')); hold on;
            H2 = PlotSwathe(100*bstSOE.bar(1:hz,ukGDP,ss), [100*bstSOE.sup2(1:hz,ukGDP,ss) 100*bstSOE.inf2(1:hz,ukGDP,ss)], rgb('dark red'));

        end
        
        handle(1)=H2.bar; handle(2) = H2.patch; handle(3)=H.patch;
        plot(zeros(nSteps),'LineWidth',0.5,'LineStyle','-','Color',rgb('black')); hold on
        
        if cc == 1
            
            plot(100*squeeze(bstUS.barCF(1:hz,usGDP,ss)),'LineWidth',2,'LineStyle','--','Color',rgb('black')); hold on
            
            
        else
            
            plot(100*bstSOE.barCF(1:hz,ukGDP,ss),'LineWidth',2,'LineStyle','--','Color',rgb('black')); hold on

        end
        
        if cc == 1
            
            title([sNames(ss) ' (US)'])
            
        else
            
            title([sNames(ss) ' (UK)'])
            
        end
            
        ylabel(lNames(usGDP))
        xlabel('Months')
        axis tight
        set(gca,'Layer','top')
        optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)
        index = index+1;
        
    end
       
end

opt=LegOption; opt.handle = handle; LegSubplot({'Mean','68% C.I.','90% C. I.'},opt);
optfont = FigFontOption(fontSize); optfont.title_weight = 'bold'; FigFont(optfont)
SaveFigure([folderOut 'counterfactual GDP'],highQuality)

close all

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
        plot(100*squeeze(bstUS.barCF(1:hz,hh,ss)),'LineWidth',2,'LineStyle','--','Color',rgb('black')); hold on
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
    SaveFigure([folderOut 'US ' sNames{ss} ' counterfactual'],highQuality)
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
        plot(100*bstSOE.barCF(1:hz,hh,ss),'LineWidth',2,'LineStyle','--','Color',rgb('black')); hold on
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
    SaveFigure([folderOut 'UK ' sNames{ss} ' counterfactual'],highQuality)
    clf('reset')

end

close all
