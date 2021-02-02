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


% minLags_x = [-1,-2, 2];
% maxLags_x = [1, 3, 5];
%
% minLags_y = [-1,-3, 1];
% maxLags_y = [1, 2, 2];

maxDev = [0,0,0];
maxDev_pred = 0;

for n = [8,10,15]
    
    % some data
    x = [(1:n)',((n+1):(2*n))'];
    %     x = [(1:n)'];
    y = (1:(3*n))';
    
    for iLag_x = 1:3
        minLag_x = minLags_x(iLag_x);
        maxLag_x = maxLags_x(iLag_x);
        
        for iLag_y = 1:3
            minLag_y = minLags_y(iLag_y);
            maxLag_y = maxLags_y(iLag_y);
            
            minLag = min(minLag_x,minLag_y);
            maxLag = max(maxLag_x,maxLag_y);
            
            for removeMean = [false,true]
                
                for padded = [true] % TODO unpad
                    
                    for iB = [1,n+1]
                        
                        %%                       
                        % [1] --- forming XtX and Xty directly
                        % size determined by x
                        
                        opt = struct();
                        opt.removeMean = removeMean;
                        opt.unpad.do = ~padded;
                        
                        dims = LM.CCAdims(size(x,1),size(y,1), iB, minLag_x, maxLag_x, minLag_y, maxLag_y);
                        
                        opt.nPad = maxLag_x - minLag_x + dims.X.nZerosTop + dims.X.nZerosBottom;
                        [XtX_FFT,xF,mX,Xtop,Xbottom,n_m] = LM.laggedXtX(x,minLag_x,maxLag_x,opt);
                        
                        opt.nPad = maxLag_y - minLag_y + dims.Y.nZerosTop + dims.Y.nZerosBottom;
                        [YtY_FFT,yF,mY,Ytop,Ybottom,n_m_] = LM.laggedXtX(y(dims.yb:dims.ye),minLag_y,maxLag_y,opt);
                        
                        assert( n_m == n_m_ );
                        
                        opt.dims = dims;
                        XtY_FFT = LM.laggedXtY(xF,yF,minLag_x,maxLag_x,minLag_y,maxLag_y,n_m,...
                            mX,Xtop,Xbottom,...
                            mY,Ytop,Ybottom,...
                            opt);
                        
                        
                        %% or
                        % [2] --- forming X and Y matrices and then XtX & Xty
                        if padded
                            [xp,yp] = LM.padCCA(x, y, iB, minLag_x, maxLag_x, minLag_y, maxLag_y);
                        else
                            xp = x;
                            yp = y;
                        end
                        
                        X = LM.lagMatrix(xp, maxLag_x - minLag_x + 1);
                        Y = LM.lagMatrix(yp, maxLag_y - minLag_y + 1);
                        
                        if removeMean
                            X = X - mean(X,1);
                            Y = Y - mean(Y,1);
                        end
                        
                        XtX = X' * X;
                        YtY = Y' * Y;
                        XtY = X'*Y;
                        
                        
                        %%
                        % checks
                        % --- sanity checks
                        assert(size(X,1) == size(Y,1));
                        if ~removeMean
                            % when removeMean the zero padding will not be
                            % visible
                            zeroRows_X = all(X == 0,2);
                            zeroRows_Y = all(Y == 0,2);
                            
                            assert( ~any( zeroRows_X & zeroRows_Y ) );
                            assert( dims.X.nZerosTop == (find(~zeroRows_X,1)-1));
                            assert( dims.X.nZerosBottom == (size(X,1) - find(~zeroRows_X,1,'last')));
                            
                            assert( dims.Y.nZerosTop == (find(~zeroRows_Y,1)-1));
                            assert( dims.Y.nZerosBottom == (size(Y,1) - find(~zeroRows_Y,1,'last')));
                        end
                        
                        assert( all(size(XtX) == size(XtX_FFT)) );
                        assert( all(size(YtY) == size(YtY_FFT)) );
                        assert( all(size(XtY) == size(XtY_FFT)) );
                        
                        
                        % --- equivalence checks
                        % these should be the same
                        m = abs((XtX-XtX_FFT) ./ XtX);
                        m = m(XtX ~= 0 ); % avoiding division by 0
                        maxDev(1) = max( maxDev(1), max(m));
                        
                        m = abs((YtY-YtY_FFT) ./ YtY);
                        m = m(YtY ~= 0 );
                        maxDev(2) = max( maxDev(2), max(m));
                        
                        m = abs((XtY-XtY_FFT) ./ XtY);
                        m = m(XtY ~= 0 );
                        maxDev(3) = max( maxDev(3), max(m));
                        
                        if maxDev(3) > 0.1
                            return;
                        end
                        
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

minLag_y = -1;
maxLag_y = 5;

x = [(1:5)'];%,(11:15)'];
y = (1:20)';

% x = (1:8)';
% y = (1:(3*24))';

%
iB = 11;

padded = true;
removeMean = false;

minLag = min(minLag_x,minLag_y);
maxLag = max(maxLag_x,maxLag_y);


if padded
    %     xp = LM.padCCAx(x, minLag_x, maxLag_x, minLag, maxLag);
    %     yp = LM.padCCAx(y, minLag_y, maxLag_y, minLag, maxLag);
    
    [xp,yp] = LM.padCCA(x, y, iB, minLag_x, maxLag_x, minLag_y, maxLag_y);
else
    
end

X = LM.lagMatrix(xp, maxLag_x - minLag_x + 1);
Y = LM.lagMatrix(yp, maxLag_y - minLag_y + 1);

%
opt = struct();

opt.removeMean = removeMean;
opt.unpad.do = ~padded;

dims = LM.CCAdims(size(x,1),size(y,1), iB, minLag_x, maxLag_x, minLag_y, maxLag_y);

opt.nPad = maxLag_x - minLag_x + dims.X.nZerosTop + dims.X.nZerosBottom;
[XtX_FFT,xF,mX,Xtop,Xbottom,n_m] = LM.laggedXtX(x,minLag_x,maxLag_x,opt);

opt.nPad = maxLag_y - minLag_y + dims.Y.nZerosTop + dims.Y.nZerosBottom;
[YtY_FFT,yF,mY,Ytop,Ybottom,n_m_] = LM.laggedXtX(y(dims.yb:dims.ye),minLag_y,maxLag_y,opt);

opt.dims = dims;
XtY_FFT = LM.laggedXtY(xF,yF,minLag_x,maxLag_x,minLag_y,maxLag_y,n_m,...
    mX,Xtop,Xbottom,...
    mY,Ytop,Ybottom,...
    opt);

XtY = X'*Y;
X,Y,XtY
int64(X'*X) - int64(XtX_FFT)
int64(Y'*Y) - int64(YtY_FFT)
int64(XtY) - int64(XtY_FFT)
% xp'
% yp'
% size(X), size(Y)
% [X,Y]
