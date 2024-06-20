function [rgb] = func_hyperImshow( hsi, RGBbands )
%% Hyperspectral Image color display
% Author: Zephyr Hou
% Time: 2019-12-02
% Function Usage
% Input:
%    hsi—the 3D hyperspectral dataset with the size of rows x cols x bands
%    RGBbands— the RGB bands to be displayed, with the format [R G B]
% Output:
%    rgb- the finual result with the RGB bands with the size of (rows x cols x 3)  

% 使用方法如下：
% load('Viareggio.mat') % 导入需要显示的影像
% [rgb] = func_hyperImshow(hsi,[219,144,66]); % func_hyperImshow(3D高光谱图像,[R,G,B]) 三个数波段号
% https://blog.csdn.net/NBDwo/article/details/115564548
%% Main Function

hsi=double(hsi);
[rows, cols, bands] = size(hsi);

minVal =min(hsi(:));
maxVal=max(hsi(:));
normalizedData=hsi-minVal;

if(maxVal==minVal)
    normalizedData=zeros(size(hsi));
else
    normalizedData=normalizedData./(maxVal-minVal);
end

hsi=normalizedData;

[rows, cols, bands] = size(hsi);

if (nargin == 1)
    RGBbands = [bands round(bands/2) 1];
end

if (bands ==1)
    red = hsi(:,:);
    green = hsi(:,:);
    blue = hsi(:,:);
else
    red = hsi(:,:,RGBbands(1));
    green = hsi(:,:,RGBbands(2));
    blue = hsi(:,:,RGBbands(3));
end

rgb = zeros(size(hsi, 1), size(hsi, 2), 3);
rgb(:,:,1) = adapthisteq(red);      % Adaptive histogram equalization
rgb(:,:,2) = adapthisteq(green);
rgb(:,:,3) = adapthisteq(blue);
imshow(rgb); axis image;   % 保持图像的显示比例，其中axis为坐标轴的控制函数。
end
