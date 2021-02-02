% Testing function part of the Linear Model (LM) package.
% Author: Octave Etard
%
% This script tests the equivalence between forming the X and y matrcies
% and then computing XtX = X' * X and Xty = X' * y, and computing XtX and
% Xty directly.
%
minLags = [-5,-2, 2];
maxLags = [-2, 3, 5];

maxDev = [0,0];
maxDev_pred = 0;

tic;
for n = [8,10,15]
    
    % some data
    x = [(1:n)',((n+1):(2*n))'];
    y = (1:(3*n))';
    
    for iLag = 1:3
        minLag = minLags(iLag);
        maxLag = maxLags(iLag);
        
        for removeMean = [true,false]
            
            for padded = [false,true]
                
                for precomputeYFFT = [false,true]
                    
                    for iB = [1,n+1]
                        
                        %%
                        % [1] --- forming XtX and Xty directly
                        opt = struct();
                        opt.removeMean = removeMean;
                        
                        if ~padded
                            % only required if unpadding
                            opt.unpad = LM.laggedDims(size(x,1),iB,size(y,1),minLag,maxLag);
                        end
                        opt.unpad.do = ~padded;
                        
                        [XtX_FFT,xF,mX,Xtop,Xbottom] = LM.laggedXtX(x,minLag,maxLag,opt);
                        opt.iB = iB;
                        opt.nx = size(x,1);
                        
                        if precomputeYFFT
                            [y_,mY,n_mY,Ytop,Ybottom] = LM.computeYFFT(y,minLag,maxLag,size(xF,1),opt);
                        else
                            y_ = y;
                            n_mY = [];
                            mY = [];
                            Ytop = [];
                            Ybottom = [];
                        end
                        
                        Xty_FFT = LM.laggedXty(xF,y_,minLag,maxLag,...
                            mX,Xtop,Xbottom,...
                            precomputeYFFT,mY,n_mY,Ytop,Ybottom,...
                            opt);
                        
                        %% or
                        % [2] --- forming X and Y matrices and then XtX & Xty
                        if padded
                            [xp,yp] = LM.pad(x,y,minLag,maxLag,iB);
                            iB_ = 1;
                        else
                            xp = x;
                            yp = y;
                            iB_ = iB;
                        end
                        
                        X = LM.laggedX(xp,minLag,maxLag,iB_,size(yp,1));
                        Y = LM.laggedY(yp,minLag,maxLag,iB_,size(xp,1));
                        
                        if removeMean
                            X = X - mean(X,1);
                            Y = Y - mean(Y,1);
                        end
                        
                        XtX = X' * X;
                        Xty = X' * Y;
                        
                        %%
                        assert( all(size(XtX_FFT) == size(XtX)) );
                        assert( all(size(Xty) == size(Xty_FFT)) );
                        
                        % these should be the same
                        maxDev(1) = max( maxDev(1), max(abs((XtX-XtX_FFT) ./ XtX),[],'all'));
                        maxDev(2) = max( maxDev(2), max(abs((Xty-Xty_FFT) ./ XtX),[],'all'));
                        
                        
                        %% testing predictions
                        nLags = maxLag - minLag + 1;
                        nFeatures = 2;
                        % some coeffs, does not matter
                        coeffs = (1:(nFeatures*nLags))';
                        
                        % through matrix multiplication
                        pred = X * coeffs;
                        
                        % through convolution
                        coeffs = reshape(coeffs,[nLags,nFeatures]);
                        if removeMean
                            m = sum(reshape(mX,nLags,nFeatures) .* coeffs,1);
                        else
                            m = 0;
                        end
                        
                        if padded
                            pred_fft = squeeze(sum(LM.convfft(x,coeffs,'full') - m,2));
                        else
                            xb = opt.unpad.xb;
                            xe = opt.unpad.xe;
                            x_ = x(xb:xe,:);
                            pred_fft = squeeze(sum(LM.convfft(x_,coeffs,'valid') - m,2));
                        end
                        maxDev_pred = max(maxDev_pred, max(abs(pred-pred_fft),[],'all'));

                    end
                end
            end
        end
    end
end
toc;
maxDev
maxDev_pred