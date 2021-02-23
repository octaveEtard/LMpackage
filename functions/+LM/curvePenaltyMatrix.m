function L = curvePenaltyMatrix(nLags,nFeatures,inverse)
%
% LM.curvePenaltyMatrix
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
%
% 2TD penalty. Actually we care mostly about inv(L) which can also be computed analytically
% and is also sparse by block.
%

% assuming X is in the form 'Feature_n x nLags ... Feature_n n nLags'
Lblock = zeros(nLags,nLags);

% --- make Lblock ---
if inverse
    % first column
    Lblock(1:(nLags-1)) = (nLags-1):-1:1;
    % last column
    Lblock(2:nLags,nLags) = 1:(nLags-1);
    
    for iStep = 1:(nLags/2)
        val = -( (iStep*iStep):iStep:((nLags-iStep-1)*iStep) );
        
        Lblock((iStep+1):(nLags-iStep),nLags-iStep) = val;
        Lblock(nLags-iStep,(iStep+1):(nLags-iStep)) = val;
        
        Lblock((nLags-iStep):-1:(iStep+1),iStep+1) = val;
        Lblock(iStep+1,(nLags-iStep):-1:(iStep+1)) = val;
    end
    
    Lblock = Lblock ./ (nLags-1);
    
else
    % diagonal
    Lblock([1,nLags*nLags]) = 1;
    Lblock((nLags+2):(nLags+1):(nLags*nLags-nLags-1)) = -2;
    % lower diagonal
    Lblock(2:(nLags+1):(nLags*nLags-2*nLags-1)) = 1;
    % upper diagonal
    Lblock((2*nLags+2):(nLags+1):(nLags*nLags-1)) = 1;
    % ----
end
% --- fill L ---
if 1 < nFeatures
    L = zeros(nFeatures * nLags,nFeatures * nLags);
    idx = 1:nLags;
    for iFeature = 1:nFeatures
        L(idx + (iFeature-1)*nLags,idx + (iFeature-1)*nLags) = Lblock;
    end
else
    L = Lblock;
end
end
%
%
