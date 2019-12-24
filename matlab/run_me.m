%%% Script to call the LC segmentation for some illustrative examples 

% ------------------------------------------------------------------------
for k=1:4

% Loading the #k-th example
load (strcat('Example_', num2str(k)));

% Computing the segmentation
[~, Ibin] = LCseg(Iorig, maskconstraints);

% Cutting the segmentation
Icut = LCcut(Iorig, Ibin, 200);

% Printing the images
disp('Printing the result');
LCoutput(Imarked, Icut);

end
% ------------------------------------------------------------------------