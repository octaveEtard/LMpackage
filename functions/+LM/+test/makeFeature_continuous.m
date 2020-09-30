function feature = makeFeature_continuous(dur,Fs,tWin)
%
% LM.test.makeFeature_continuous
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
nPnts = ceil(dur*Fs)+1;
feature = abs(randn(nPnts,1));
feature = feature / max(feature);

if 0 < tWin
    nWin = ceil(tWin * Fs) + 1;
    win = hann(2*nWin-1);
    feature(1:nWin) = win(1:nWin) .* feature(1:nWin);
    feature((1:nWin)-nWin+end) = win(nWin:end) .* feature((1:nWin)-nWin+end);
end

end
%
%