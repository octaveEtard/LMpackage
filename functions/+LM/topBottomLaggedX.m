function [Xtop,Xbottom] = topBottomLaggedX(x,xb,xe,minLag,maxLag)
%
% topBottomLaggedX
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%

[nPnts,nFeatures] = size(x);
% [xb,xe] = LM.laggedDimsX(nPnts,minLag,maxLag);

nLags = maxLag - minLag + 1;
nTop = xb - 2 + nLags;
nBottom = nLags + nPnts - xe - 1;


%%
% Xtop = nan(nTop, nLags * nFeatures,'like',x);
% Xbottom = nan(nBottom, nLags * nFeatures,'like',x);
% 
% idxTop = (1:nTop) + xb - 1;
% idxBottom = (1:nBottom) + xe - nBottom - nLags + 1;
% 
% dimOffset = ((1:nFeatures)-1)*nLags; 
% 
% for iLag = 1:nLags
%     Xtop(:,dimOffset + iLag) = x(idxTop - iLag + nLags,:);
%     Xbottom(:,dimOffset + iLag) = x(idxBottom - iLag + nLags,:);
% end


%%
Xtop = zeros(nTop, nLags * nFeatures,'like',x);
Xbottom = zeros(nBottom, nLags * nFeatures,'like',x);

dimOffset = ((1:nFeatures)-1)*nLags; 

for iLag = 1:nLags
    
    n = nTop - iLag + 1;
    Xtop((1:n) + iLag - 1,dimOffset + iLag) = x(1:n,:);
    
    n =  nPnts - xe + iLag - 1 ;
    Xbottom(1:n,dimOffset + iLag) = x((1:n)-n+end,:);
end



end