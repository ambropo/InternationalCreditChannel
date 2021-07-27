function out = cookVARData(DATA,obs,nvar,fo,lo,XoXlag,vTreat,vNames)
out = nan(obs,nvar);
for ii=1:nvar
    if vTreat(ii)==0
        out(:,ii) = DATA.(vNames{ii})(fo:lo,1);
    elseif vTreat(ii)==1
        out(:,ii) = log(DATA.(vNames{ii})(fo:lo,ii));
    elseif vtreat(ii)==2
        out(:,ii) = XoX(DATA.(vNames{ii})(fo:lo,1),XoXlag,'logdiff');
    elseif vtreat(ii)==3
        out(:,ii) = XoX(DATA.(vNames{ii})(fo:lo,1),XoXlag,'diff');
    end
end

end