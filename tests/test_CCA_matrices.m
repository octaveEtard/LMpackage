% Testing function part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Testing function part of the Linear Model (LM) package.
% Author: Octave Etard
%
% This script tests the equivalence between forming the X and y matrcies
% and then computing XtX = X' * X and Xty = X' * y, and computing XtX and
% Xty directly.
%
minLags_x = [-5,-2, 2];
maxLags_x = [-2, 3, 5];

minLags_y = [-6,-3, 1];
maxLags_y = [-1, 2, 2];

maxDev = [0,0];
maxDev_pred = 0;

for n = [8,10,15]
    
    % some data
    x = [(1:n)',((n+1):(2*n))'];
    y = (1:n)';
    
    for iLag_x = 1:3
        minLag_x = minLags_x(iLag_x);
        maxLag_x = maxLags_x(iLag_x);
        
        for iLag_y = 1:3
            minLag_y = minLags_y(iLag_y);
            maxLag_y = maxLags_y(iLag_y);
            
            minLag = min(minLag_x,minLag_y);
            maxLag = max(maxLag_x,maxLag_y);
            
            for removeMean = [false]
                
                for padded = [true]
                    
                    for iB = [1,n+1]
                        
                        %%
                        % [1] --- forming XtX and Xty directly
                        opt = struct();
                        opt.nPad = maxLag - minLag;
                        opt.removeMean = removeMean;
                        opt.unpad.do = ~padded;
                        
                        XtX_FFT = LM.laggedXtX(x,minLag_x,maxLag_x,opt);
                        YtY_FFT = LM.laggedXtX(y,minLag_y,maxLag_y,opt);
                        
                        %% or
                        % [2] --- forming X and Y matrices and then XtX & Xty
                        if padded
                            xp = LM.padCCA(x, minLag_x, maxLag_x, minLag, maxLag);
                            yp = LM.padCCA(y, minLag_y, maxLag_y, minLag, maxLag);
                        else
                            error('!');
                        end
                        
                        X = LM.lagMatrix(xp, maxLag_x - minLag_x + 1);
                        Y = LM.lagMatrix(yp, maxLag_y - minLag_y + 1);
                        
                        if removeMean
                            error('!');
                        end
                        
                        XtX = X' * X;
                        YtY = Y' * Y;
                        
                        %%
                        assert( all(size(XtX) == size(XtX_FFT)) );
                        assert( all(size(YtY) == size(YtY_FFT)) );
                        
                        % these should be the same
                        maxDev(1) = max( maxDev(1), max(abs((XtX-XtX_FFT) ./ XtX),[],'all'));
                        maxDev(2) = max( maxDev(2), max(abs((YtY-YtY_FFT) ./ YtY),[],'all'));
                        
                        
                        %% testing predictions
                        
                    end
                end
            end
        end
    end
end
maxDev
maxDev_pred



return
%%
minLag_x = -1;
maxLag_x = 1;

minLag_y = -2;
maxLag_y = 2;


minLag = min(minLag_x,minLag_y);
maxLag = max(maxLag_x,maxLag_y);

x = (1:5)';
y = (6:10)';

xp = LM.padCCA(x, minLag_x, maxLag_x, minLag, maxLag);
yp = LM.padCCA(y, minLag_y, maxLag_y, minLag, maxLag);

X = LM.lagMatrix(xp, maxLag_x - minLag_x + 1);
Y = LM.lagMatrix(yp, maxLag_y - minLag_y + 1);

X,Y
X'*Y

% xp'
% yp'
% size(X), size(Y)
% [X,Y]

%%
opt = [];
opt.nPad = maxLag - minLag;
opt.unpad.do = false;
opt.removeMean = false;

opt.minLag_x = minLag_x;
opt.maxLag_x = maxLag_x;

opt.minLag_y = minLag_y;
opt.maxLag_y = maxLag_y;

[XtX,xF] = LM.laggedXtX(x,minLag_x,maxLag_x,opt);
[yty,yF] = LM.laggedXtX(y,minLag_y,maxLag_y,opt);

[XtY] = LM.laggedXt_2(xF,yF,opt)