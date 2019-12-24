function main_interactive(imgName)
% ------------------------------------------------------------------------
% Description: interactive main function (example of user interface)
% Usage: main_interactive('sample_image.png');
% 
% Input(s):
%   imgName: Filename of the target image + extension
%
% Output(s):
%   Segmentation computed by the LC framework
% ------------------------------------------------------------------------

% Loading image/var initialization
Iorig = imread(imgName);
Imarked = Iorig;
[m, n, ~] = size(Iorig);
maskconstraints = zeros(m,n,2);
colorconstraints = zeros(2,3);
first = true;

while (true)

  % Generating the mask/brushed image
  [mask, Imarked, color] = LCseeds(Imarked);
   
  % Checking the seed type (notice that only foreground and background 
  % seeds are supported by this script)
  if ((sum(color(1,:) == colorconstraints(1,:)) == 3) || (first))
      k = 1; l = 2;
  else
      k = 2; l = 1;
  end
  
  maskconstraints(:,:,k) = logical(maskconstraints(:,:,k) + mask(:,:,1));
  maskconstraints(:,:,l) = logical(maskconstraints(:,:,l) + mask(:,:,2));    
  colorconstraints = color;
  first = false;
  
  % Computing the LC-based segmentation
  [~, Ibin] =  LCseg(Iorig, maskconstraints);

  % Cutting the segmentation
  Icut = LCcut(Iorig, Ibin, 200);

  % Printing the images
  disp('Printing the result');
  fig = LCoutput(Imarked, Icut);
  
  disp('Press any key to continue or click any mouse button to quit');
  button = waitforbuttonpress;
  close(fig);
  
  % If the mouse is activated, then close the application 
  if (button == 0)
    imwrite(mat2gray(Imarked), strcat(imgName(1:end-4),'_marked.png'));
    imwrite(mat2gray(Ibin), strcat(imgName(1:end-4),'_part.png'));
    imwrite(mat2gray(Icut), strcat(imgName(1:end-4),'_cut.png'));
    imwrite(mat2gray(maskconstraints(:,:,1)), strcat(imgName(1:end-4),'_mask1.png'));
    imwrite(mat2gray(maskconstraints(:,:,2)), strcat(imgName(1:end-4),'_mask2.png'));
    break;
  end
  
end

end