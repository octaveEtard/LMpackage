function [impResponse,t,nMax] = makeIR(tMax,Fs,amplitude,delay,sigma)
%
% LM.test.makeIR
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Make synthetic impulse response as a sum of gaussian pdfs
%

nMax = ceil(tMax*Fs)+1;
t = (-nMax:nMax)'/Fs;
nPulses = numel(amplitude);

impResponse = zeros(2*nMax+1,1);

for iPulse = 1:nPulses
    impResponse = impResponse + amplitude(iPulse) * ...
        sigma(iPulse) * sqrt(2*pi) * normpdf(t,delay(iPulse),sigma(iPulse));
end

end
%
%