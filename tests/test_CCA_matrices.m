% Testing function part of the Linear Model (LM) package.
% Author: Octave Etard
%

minLag_x = -5;
maxLag_x = 1;

minLag_y = -2;
maxLag_y = 2;


minLag = min(minLag_x,minLag_y);
maxLag = max(maxLag_x,maxLag_y);

x = (1:10)';
y = (1:10)';

xp = LM.padCCA(x, minLag_x, maxLag_x, minLag, maxLag);
yp = LM.padCCA(y, minLag_y, maxLag_y, minLag, maxLag);

X = LM.lagMatrix(xp, maxLag_x - minLag_x + 1);
Y = LM.lagMatrix(yp, maxLag_y - minLag_y + 1);

xp'
yp'
size(X), size(Y)
[X,Y]
