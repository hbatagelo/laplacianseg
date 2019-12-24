function Icut = LCcut(Iorig, Ibin, bgcolor)

% Getting the image size
[m, n, c] = size(Iorig);

% Cutting the segmentation
Icut = bgcolor*ones(m,n,c);
for i=1:c
  aux = Iorig(:,:,i);
  aux(~Ibin) = bgcolor;
  Icut(:,:,i) = aux;
end

end