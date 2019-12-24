function fig = LCoutput(Imarked, Icut)

fig = figure;
subplot(1,2,1), subimage(mat2gray(Imarked));
title(strcat('Seeded image'));
axis off;

subplot(1,2,2), subimage(mat2gray(Icut));
title(strcat('Segmentation'));
axis off;

end