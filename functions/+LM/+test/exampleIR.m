function [impResponse,t,nMax] = exampleIR(type,Fs)
%
% LM.test.exampleIR
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Two examples of impulse responses.
%

tMaxResponse = 400e-3; % in s

switch type
    case 1  % simple Gaussian response
        delayResponse = 200e-3; % in s
        widthResponse = 50e-3; % in s
        amplitudeResponse = 1;
        
    case 2  % bi-modal response:
        amplitudeResponse = [1,-1.5] / 1.5;
        delayResponse =  [80,180]*1e-3; % in s
        widthResponse = [20,40]* 1e-3; % in s
        
    otherwise
        error('Condition not implemented');
end

[impResponse,t,nMax] = LM.test.makeIR(tMaxResponse,Fs,amplitudeResponse,delayResponse,widthResponse);

end
%
%