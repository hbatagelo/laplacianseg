function Feats =  LCfeats(Iorig)
% ------------------------------------------------------------------------
% Description: computes image features
% Usage: Feats = LCfeats(Iorig) 
% 
% Input(s):
%   Iorig: target image
% Output(s):
%   Feats: image features
% ------------------------------------------------------------------------

% Initializations
mm = ones(3)/(3)^2;
rs = 5;

% Preprocessing
%op : 1
%for i=1:size(Iorig,3)
%   u = imadjust(LCdiff(double(Iorig(:,:,i)), 22, 0.17, 10, 1));
%   Iorig(:,:,i) = imadjust(u);
%end

%op : 2
for i=1:size(Iorig,3)
   u = conv2(Iorig(:,:,i), mm, 'same');
   Iorig(:,:,i) = u;     
end

% Standard color features
Feat1 = uint8(applycform(mat2gray(Iorig), makecform('srgb2lab')));
Feat2 = Iorig;

SA = rangefilt(Feat1, true(rs));
Feats(:,:,1) = imadjust(uint8(SA(:,:,1)));
Feats(:,:,2) = imadjust(uint8(SA(:,:,2)));
Feats(:,:,3) = imadjust(uint8(SA(:,:,3)));

SB = rangefilt(Feat2, true(rs));
Feats(:,:,1+end) = imadjust(uint8(SB(:,:,1)));
Feats(:,:,1+end) = imadjust(uint8(SB(:,:,2)));
Feats(:,:,1+end) = imadjust(uint8(SB(:,:,3)));

end