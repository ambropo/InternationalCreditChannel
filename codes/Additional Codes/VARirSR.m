function [IRFsr,MT,p05,p95,med,MTidx] = VARirSR(VAR,VARopt,s,sigmaHat,signRestMat,nPeriods,maxIters)

nVar = VAR.nvar;
nSteps = VARopt.nsteps;
signSize = size(signRestMat,1)/(nPeriods*nVar);
sign0 = signRestMat(1:nVar,:);
warning off all
irf0 = signRestrictions(sigmaHat,s,sign0,maxIters);
nIrf0 = size(irf0,1)/nVar;
mpsr = zeros(nIrf0,1);
IRFsr = cell(nIrf0,1);
if nIrf0>0
    
    for ii = 1:nIrf0
        
        irf0temp = irf0((ii-1)*nVar+1:ii*nVar,:);
        irfTemp = VARirfivSR(VAR,VARopt,irf0temp);
        mpsr(ii) = 1;
        
        if nPeriods > 1
            
            for pp = 2:nPeriods
                
                if signSize==1
                    
                    signMat = signRestMat((pp-1)*nVar+1:pp*nVar,:);
                    
                else
                    
                    signMat = signRestMat;
                    
                end
                
                checkSign = squeeze(irfTemp(pp,:,size(s,2)+1:nVar));
                
                if sign(checkSign(logical(abs(signMat)))) == ...
                        signMat(logical(abs(signMat)))
                    
                    mpsr(ii) = mpsr(ii)+1;
                    
                else
                    
                    break
                    
                end
                
            end
            
        end
        
        if mpsr(ii) == nPeriods
            
              IRFsr{ii} = irfTemp;
            
        end
        
    end
    
else
    
    MTidx = 0;
    MT = [];
    p05 = [];
    p95 = [];
    med = [];
    
end

if sum(mpsr < nPeriods)==nIrf0
    
    MTidx = 0;
    MT = [];
    p05 = [];
    p95 = [];
    med = [];
    
else
    
        IRFsr = permute(cat(4,IRFsr{mpsr == nPeriods}),[4 1 2 3]);
        nModels = sum(mpsr == nPeriods);
%     disp(['Bingo! Found ' num2str(nModels) ' rotations for this draw']);
    IRFsrmt = IRFsr - repmat(median(IRFsr,1),[nModels 1 1 1]);
    IRFsrmtvec = zeros(nModels,nSteps);
    
    for ii = 1:nModels
        
        for jj = 1:nSteps
            
            IRFsrmtvec(ii,jj) = sum(sum(squeeze(IRFsrmt(ii,jj,:,:)).^2));
            
        end
        
    end
    
    [~,MTidx] = min(sum(IRFsrmtvec,2));
    MT = squeeze(IRFsr(MTidx,:,:,:));
    p05 = squeeze(quantile(IRFsr,.05,1));
    p95 = squeeze(quantile(IRFsr,.95,1));
    med = squeeze(quantile(IRFsr,.5,1));
    
end

end
