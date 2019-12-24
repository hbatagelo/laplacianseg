function W = LCweights(Iorig)
% ------------------------------------------------------------------------
% Description: builds the matrix of weights
% Usage: W = LCweights(I);
% 
% Input(s):
%   Iorig: original image
%
% Output(s):
%   W: Matrix of weights (sparse/adjacency matrix)
% ------------------------------------------------------------------------

% Computing the image features
I = LCfeats(Iorig);

%%% Initializations
Epsilon = 0.0001;  % range: [10^-6,10^-1]
Beta = 600; % range: [50,700]
[nRows, nCols, nChannels] = size(I);
Img = zeros(nRows, nCols, nChannels);
Dist_lr_aux = zeros(nRows-1, nCols-1, nChannels);
Dist_rl_aux = zeros(nRows-1, nCols-1, nChannels);

%%% Computing the weight intensities of the edges
% Computing the diagonal elements
for k=1:nChannels
   Img(:,:,k) = mat2gray(I(:,:,k));
   Dist_lr_aux(:,:,k) = abs(Img(1:end-1,1:end-1,k) - Img(2:end,2:end,k));
   Dist_rl_aux(:,:,k) = abs(Img(1:end-1,2:end,k) - Img(2:end,1:end-1,k));
end

% Computing the max for the image feature differences
dist_lr = reshape(max(Dist_lr_aux, [], 3), (nRows-1)*(nCols-1), 1);
dist_rl = reshape(max(Dist_rl_aux, [], 3), (nRows-1)*(nCols-1), 1);
dist_ht = reshape(max(abs(diff(Img,1,1)), [], 3), (nRows-1)*nCols, 1);
dist_vt = reshape(max(abs(diff(Img,1,2)), [], 3), nRows*(nCols-1), 1);

% Computing the weights
dist = [dist_lr; dist_rl; dist_ht; dist_vt];
dist = roundn(exp(-Beta*dist.*dist),-4) + Epsilon;

%%% Building the indices of the above weights
iaux = transpose(1:nRows*nCols);

%"Left-to-right" edge indices
ilr = iaux(1:end-nRows);
ilr(nRows:nRows:end) = [];
jlr = ilr + nRows + 1;

%"Right-to-left" edge indices
irl = ilr + 1;
jrl = jlr - 1;

%"Horizontal" edge indices
iht = iaux;
iht(nRows:nRows:end) = [];
jht = iht + 1;

%"Vertical" edge indices
ivt = iaux(1:end-nRows);
jvt = iaux(nRows+1:end);

%Setting the i and j-index vectors
i = [ilr; irl; iht; ivt];
j = [jlr; jrl; jht; jvt];

%%% Building the matrix of weights (W)
W = sparse([i; j], [j; i], [dist; dist], nRows*nCols, nRows*nCols, 2*numel(i));

end