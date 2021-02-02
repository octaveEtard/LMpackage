function opt = CCAdims(nx, ny, iB, minLag_x, maxLag_x, minLag_y, maxLag_y)
%
% CCAdims
% Part of the Linear Model (LM) package.
% Author: Octave Etard
%
% Compute useful dimensions / indices for the CCA matrices, in particular
% with respect ot zero-padding.
%
% opt.(xb/yb) : index of first sample of timeserie x/y to use
% opt.(xe/ye) : index of last sample of timeserie x/y to use
% opt.(x/y).nPad_b : number of zeros to add at the begining of the timeseries
% opt.(x/y).nPad_e : number of zeros to add at the end of the timeseries
% opt.(X/Y).nZerosTop : number of zero rows on top of resulting matrices
% opt.(X/Y).nZerosBottom : number of zero rows on bottom of resulting matrices
% opt.y.iB : index of the sample corresponding to sample 1 of x
%
minLag = min(minLag_x,minLag_y);
maxLag = max(maxLag_x,maxLag_y);

opt = [];

% --- for x
% timeserie:
opt.xb = 1;
opt.xe = nx;
opt.x.nPad_b = maxLag - minLag_x;
opt.x.nPad_e = maxLag_x - minLag;
% X matrix
opt.X.nZerosTop = maxLag - maxLag_x;
opt.X.nZerosBottom = minLag_x- minLag;

% --- for y
nPad = maxLag_y - minLag_y;

yb = iB - maxLag + maxLag_y;
if yb < 1
    opt.y.nPad_b = nPad + 1 - yb;
    opt.Y.nZerosTop = 1 - yb;
    yb = 1;
else
    opt.y.nPad_b = nPad;
    opt.Y.nZerosTop = 0;
end

ye = iB + nx - 1 - minLag + minLag_y;
if ny < ye
    opt.y.nPad_e = nPad + ye - ny;
    opt.Y.nZerosBottom = ye - ny;
    ye = ny;
else
    opt.y.nPad_e = nPad;
    opt.Y.nZerosBottom = 0;
end

opt.yb = yb;
opt.ye = ye;

opt.y.iB = iB;

end
%
%