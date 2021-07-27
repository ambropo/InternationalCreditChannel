function irf0mat= signRestrictions(sigmaHat,s,...
    signRestMat,maxRot)
warning off all
% initialisation
dy = size(signRestMat,1);
ds = size(signRestMat,2);
whereToStart = 1+dy-ds;
if ~(size(s,2)==whereToStart-1)
    error('s and signRestMat must have coherent sizes')
end
nanMat = NaN*ones(dy,1);
orderIndices = 1:dy;

if whereToStart>1 
    C = chol(sigmaHat,'lower');
    q = C\s;
    for ii = whereToStart:dy
        r = randn(dy,1);
        q = [q (eye(dy)-q*q')*r/norm((eye(dy)-q*q')*r)];
    end
    startingMat = C*q;
else
    startingMat = chol(sigmaHat,'lower');
end
irf0mat = [];
counter = 0;
while counter<maxRot
    counter = counter+1;
    termaa = startingMat;
    TermA = 0;
    rotMat = eye(dy);
    rotMat(whereToStart:end,whereToStart:end) = getqr(randn(length(whereToStart:dy)));
    %rotate columns of Sq with Hausholder matrix
    terma = termaa*rotMat;
    termaa = terma;
    %check sign restrictions
    for ii = 1 : ds
        for jj = whereToStart : dy
            if isfinite(terma(1,jj))
                if sum(termaa(:,jj) .* signRestMat(:,ii) < 0) == 0
                    TermA = TermA + 1;
                    orderIndices(whereToStart-1+ii) = jj;
                    terma(:,jj) = nanMat;
                    break
                elseif sum(-termaa(:,jj) .* signRestMat(:,ii) < 0) == 0
                    TermA = TermA + 1;
                    terma(:,jj) = nanMat;
                    termaa(:,jj) = -termaa(:,jj);
                    orderIndices(ii+whereToStart-1) = jj;
                    break
                    
                end
            end
        end
    end
    if isequal(TermA,ds)
        irf0=termaa(:,orderIndices);
        irf0mat = [irf0mat;irf0];
    end
end
% if isequal(TermA,ds)
%     irf0mat=termaa(:,orderIndices);
% else
%     irf0mat = [];
% %     disp('No solution found')
% end
end


function out=getqr(a)

[q,r]=qr(a);
for i=1:size(q,1)
    if r(i,i)<0
        q(:,i)=-q(:,i);
    end
end

out=q;

end
