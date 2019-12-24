function [maskmat, Imarked, rgbbrushcolor] = LCseeds(I)
% ------------------------------------------------------------------------
% Description: produces the seeds from a color palette
% Usage: I = imread('image1.png');
%        [maskmat Imarked] = LCseeds(I);
% 
% Input(s):
%   imgFile: Target image
%
% Output(s):
%   maskmat: Matrix of seeds (mask)
%   Imarked: Brushed image
%   rgbbrushcolor: Matrix containing the selected rgb colors
% ------------------------------------------------------------------------


%%% Convertion/var initialization
[nRows, nCols, nChannels] = size(I);

Imarked = I;
if (nChannels == 1)
   for k=2:3
      Imarked(:,:,k) = I;
   end
end

first = true;
lineColor = [0, 1, 0];
lineDensity = 3;

%%% Creating toolbar/printing the options
sProg = 'Seeded Image Segmentation: Laplacian Coordinates - Insert seeds';
fig = figure('Name', sProg, 'MenuBar', 'none', 'ToolBar', 'none');

sPress = 'Press to select a seed color';
toolbar = uitoolbar(fig);
colorButton = uipushtool(toolbar, 'CData', setColor(lineColor), ...
'TooltipString', sPress, 'ClickedCallback', @setColorCallback);

imshow(Imarked);
fprintf('Click on image to mark the seeds\n');
fprintf('c: Set the selected seeds, Esc: Run the segmentation \n');

%%% Seeding process
hold on;
x = [];
y = [];
maskmat = zeros(nRows, nCols, 1);
rgbbrushcolor = zeros(1, 3);

while (true)
  
  % Identify points from the cursor with the mouse  
  [px, py, button] = ginput(1);
  
  % Getting the color from the palette
  ico = get(colorButton, 'CData');
  
  % If press 'c'
  if (button == 99)
      
     n = length(x);
     if (x == 0) 
        continue;
     end
    
     % Creating a color line based on the input points
     mask = 0;
     for i=2:n
       xa = x(i - 1);
       ya = y(i - 1);
       xb = x(i);
       yb = y(i);
       mask = mask + drawLine([xa ya], [xb yb], lineDensity, nRows, nCols);
     end
     
     % Drawing the color line on the image
     mask = (mask > 0);
     for i=1:3
        %mark = mask * ico(1,1,i) + (1 - mask);
        %Imarked(:, :, i) = uint8(mark .* double(Imarked(:, :, i)));
        mark = mask * floor(255*ico(1,1,i));
        Imarked(:, :, i) = uint8((1 - mask).* double(Imarked(:, :, i)) + mark);
     end
     
     % Updating the image
     figure(fig);
     hold off;
     imshow(Imarked);
     x = [];
     y = [];
     hold on;
     
    % Adding current mask + chosen color into the matrix of seeds;
    lineNewColor = [ico(1,1,1), ico(1,1,2), ico(1,1,3)];
    if ((sum(lineNewColor == lineColor) == 3) || (first)) % the same type of seed 
       maskmat(:,:,end) = maskmat(:,:,end) + mask;
       rgbbrushcolor(end,:) = lineNewColor;
       first = false;
    else               % if the seed type has been changed
       maskmat(:,:,end+1) = mask;
       rgbbrushcolor(end+1,:) = lineNewColor;
    end
    lineColor = lineNewColor;
 
  % If press 'Esc'  
  elseif (button == 27)
    
    close(fig);
    break;
  
  else
      
    figure(fig);
    plot(px, py, 'x');
    x = [x; px];
    y = [y; py];
  
  end
  
end

end


%%% Support functions
%%% =======================================================================
function mask = drawLine(xa , xb, lineDensity, nRows, nCols)

[xG, yG] = meshgrid(1:nCols,1:nRows);
mask = zeros(nRows,nCols);
d = xb - xa;
step = (lineDensity / norm(d));

for t=0:step:1
  xn = xa + t * d;
  dImage = (xG - xn(1)).^2 + (yG - xn(2)).^2;
  mask = mask + (dImage < lineDensity^2);
end

mask = (mask > 0);

end


function ico = setColor(color)

ico = zeros(16, 16, 3);
for i=1:3
  ico(:,:,i) = color(i);
end

end


function setColorCallback(hObject, ~)

color = uisetcolor(hObject, 'Choose the seed color');
if ~isequal(color, 0)
  set(hObject, 'CData', setColor(color));
end

end
%%% =======================================================================