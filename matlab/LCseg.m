function [Isol, Ibin] =  LCseg(Iorig, maskconstraints)
% ------------------------------------------------------------------------
% Description: computes the foreground/background LC-based segmentation
% Usage: [Iso, Ibin] = LCSeg(Iorig, maskconstraints) 
% 
% Input(s):
%   Iorig: Target image (loaded with 'imread' function)
%   maskconstraints: Foreground/Background matrix of masks (mxnx2)
% Output(s):
%   Isol: Scalar map (obtained from the linear system of equations)
%   Ibin: Foreground/Background binary segmentation
% ------------------------------------------------------------------------

% Computing the matrix of weights
disp('Computing the matrix of weights');
W = LCweights(Iorig);

% Building/solving the linear system of equations
disp('Building/solving the LC linear system');
[A, b, labels, cindex, xsol] = LClinear_hard(W, maskconstraints); 
xunk = A\b;

% Outlier pruning
disp('Segmenting the image');
xunk(xunk < labels(1)) = labels(1);
xunk(xunk > labels(2)) = labels(2);

% Assigning the obtained solution to the output vector
xsol(cindex) = xunk;
Isol = reshape(xsol, size(Iorig,1), size(Iorig,2));

% Generating the definitive segmentation
%op: 1 (by Otsu's thresholding)
Ibin = (Isol > graythresh(Isol));
%op: 2 (by trivial cutting)
%Ibin = (Isol > sum(labels)/2);

end