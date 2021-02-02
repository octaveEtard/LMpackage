function [yF,mY,n_mY,Ytop,Ybottom] = computeYFFT(y,minLag,maxLag,nFFT,opt)
%
% LM.computeYFFT
%
% Compute the conjugated of the Fourier Transform of y that is necessary to
% make X' * y.
%
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
iB = opt.iB;
nx = opt.nx;

ny = size(y,1);

yb = max(1,iB-maxLag);
ye = min(ny,iB+nx-1-minLag);

yF = conj(fft(y(yb:ye,:),nFFT,1));

n_mY = 0;
mY = [];
Ytop = [];
Ybottom = [];

if opt.unpad.do
    [Ytop,Ybottom] =  LM.topBottomLaggedY(...
        y,opt.unpad.yb_pad_top,opt.unpad.ye_pad_top,...
        opt.unpad.yb_pad_bottom,opt.unpad.ye_pad_bottom,...
        opt.unpad.nZeros_top_y,opt.unpad.nZeros_bottom_y);
end

if opt.removeMean || 1 < nargout 
    if opt.unpad.do
        n_mY = opt.unpad.ye - opt.unpad.yb + 1;
        mY = mean(y(opt.unpad.yb:opt.unpad.ye,:),1);
    else
        nLags = maxLag - minLag + 1;
        n_mY = (nx + nLags - 1);
        mY = yF(1,:) / n_mY;
    end
end


end