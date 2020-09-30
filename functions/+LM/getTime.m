function [t,coeffs] = getTime(opt,Fs,type,coeffs,timeDim)
%
% getTime
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Return a time vector, and coeff such that the causal part of the
% coefficients are at positive latencies. timeDim is the dimension
% corresponding to time in coeffs.
%
if nargin < 5
    % time dimension = first one by default
    timeDim = 1;
end
if nargin < 4
    coeffs = [];
end

switch type
    case 'forward'
        t =  (-opt.maxLag:-opt.minLag) / Fs;
    case 'backward'
        t =  (opt.minLag:opt.maxLag) / Fs;
        
        if ~isempty(coeffs)
            coeffs = flip(coeffs,timeDim);
        end
end

end
%
%