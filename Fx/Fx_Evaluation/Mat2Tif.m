% https://blog.csdn.net/jing_zhong/article/details/111770040
% 
%InputMatFileName如D:/360极速浏览器下载/Pavia.mat
%OutputTifFilename如D:/360极速浏览器下载/PaviaTIF.tif
% Mat2Tif('D:/360极速浏览器下载/Pavia.mat','D:/360极速浏览器下载/TifPavia.tif');
function Mat2Tif(InputMatFileName,OutputTifFilename)
    load(InputMatFileName);
    InputMatImg = double(output);
    t = Tiff(OutputTifFilename,'w');
    if size(InputMatImg,3) == 3
        t.setTag('Photometric',Tiff.Photometric.RGB);
    else
        t.setTag('Photometric',Tiff.Photometric.MinIsBlack);%颜色空间解释方式
    end
    t.setTag('Compression',Tiff.Compression.None);%无压缩
    t.setTag('BitsPerSample',64);% 由于输入.mat为double类型，所以选择了64位
    t.setTag('SamplesPerPixel',size(InputMatImg,3));% 每个像素的波段数目
    t.setTag('SampleFormat',Tiff.SampleFormat.IEEEFP);% 配合BitsPerSample64位double类型，选择IEEEFP来对应
    t.setTag('ImageLength',size(InputMatImg,1));% 影像宽度
    t.setTag('ImageWidth',size(InputMatImg,2));% 影像高度
    t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);%平面配置选择集中式
    t.write(InputMatImg);% 准备好了头文件，开始写影像数据
    t.close();% 关闭影像
%下面的代码仅为了测试显示结果
%tdc = Tiff('D:/360极速浏览器下载/TifPavia.tif','r');
%TifPavia = read(tdc);
%close(tdc);
% i=12;j=55;k=89;%自选三个波段
% InputImg_r= TifPavia(:,:,i);%获取第i个波段
% InputImg_g= TifPavia(:,:,j);%获取第j个波段
% InputImg_b= TifPavia(:,:,k);%获取第k个波段
% InputImg_r= uint8(InputImg_r);%将i波段的灰度值转为0~255
% InputImg_g= uint8(InputImg_g);%将j波段的灰度值转为0~255
% InputImg_b= uint8(InputImg_b);%将k波段的灰度值转为0~255
% RGBImg=cat(3,InputImg_r,InputImg_g,InputImg_b);%将i、j、k三个波段进行合成
% figure;
% subplot(221);imshow(InputImg_r);title('红色波段');
% subplot(222);imshow(InputImg_g);title('绿色波段');
% subplot(223);imshow(InputImg_b);title('蓝色波段');
% subplot(224);imshow(RGBImg);title('合成波段');
