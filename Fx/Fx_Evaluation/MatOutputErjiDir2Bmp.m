function [] = MatOutputErjiDir2Bmp (FusionImgYijiPath,Step,saveDir)
addpath(genpath('.\Toolbox\'));
    
    % 遍历出二级目录名
    ErjiDir_list = dir(FusionImgYijiPath) ;  % 二级目录列表
    ErjiDir_list_Nums = size(ErjiDir_list,1);  % 二级目录个数 包括 .和..
    for i_ErjiDir = 3 : ErjiDir_list_Nums
        %列出当前二级文件夹内所有的文件
        ErjiPath = fullfile(FusionImgYijiPath,ErjiDir_list(i_ErjiDir).name); 
        FusionImg_list = dir([ErjiPath,'\','*.mat']) ;
        % % 获取文件名并排序
        % MatNames = {FusionImg_list.name};
        % SortMatNames = sort_nat(MatNames); 
        % FusionImg_list(1,:) = SortMatNames;
        % 在当前二级目录处理每一个图像         
        NumImgs = size(FusionImg_list,1);  % mat个数 . ..
        for i = 1:NumImgs
            % 判断是否每若干张一处理的Step
            if mod(i, Step) == 0
                formatSpec = '正在处理目录 %s！%d个图像中每隔 %d 个一处理，这是第%d个！\n';
                fprintf(formatSpec,ErjiPath, NumImgs,Step, i);
                
                    
                % 然后再正常运行
                load([FusionImg_list(i).folder,'\',FusionImg_list(i).name]); 
    %             SaveName = fullfile(saveDir,[num2str(i),'.bmp']);
    
                suffix = find('.'== FusionImg_list(i).name); %寻找后缀名前面的标志"."
                imname = FusionImg_list(i).name(1:suffix-1); %取文件名第一位到.的前一位字符
                saveErjiDir = fullfile(saveDir,ErjiDir_list(i_ErjiDir).name);
                saveName = fullfile(saveErjiDir,[imname,'.bmp']);
                
                if ~exist(saveErjiDir,'dir')%待保存的图像文件夹不存在，就建文件夹
                    mkdir(saveErjiDir)            
                end
    
                % imwrite(mat2gray(output(:,:,4:-1:2)),saveName);
               
                [rgb] = func_hyperImshow(output,[1,2,3]); % func_hyperImshow(3D高光谱图像,[R,G,B])
                imwrite(rgb,saveName);

    
%                 % 进行2%的线性拉伸https://www.zhihu.com/question/526858450        
%                 % 读入灰度图像
%                 % img = imread(saveName);
%                 % img = im2double(img);%将uint8类型的数据转换为double类型的同时，把数据范围由原来的0~255映射到0~1，可以看作数据的一种归一化,以便计算
% 
%                 img = mat2gray(output(:,:,4:-1:2));
% 
% % 计算最小值和最大值
% low_RGB_triplet(1) = prctile(img(:,1), 2);
% low_RGB_triplet(2) = prctile(img(:,2), 2);
% low_RGB_triplet(3) = prctile(img(:,3), 2);
% high_RGB_triplet(1) = prctile(img(:,1), 98);
% high_RGB_triplet(2) = prctile(img(:,2), 98);
% high_RGB_triplet(3) = prctile(img(:,3), 98);
% 
%                 % 确保 low_RGB_triplet 中的元素小于 high_RGB_triplet 中的元素
%                 % 使用了min函数来确保low_RGB_triplet中的每个元素都不大于相应的high_RGB_triplet中的元素。eps是一个很小的正数，以确保两者之间有足够的差异。这样，你就可以避免imadjust函数中的错误。IMADJUST: LOW_IN 必须小于 HIGH_IN
%                 low_RGB_triplet = min(low_RGB_triplet, high_RGB_triplet - eps);
% 
%                 % 对图像进行线性拉伸
%                 img_normalized = mat2gray(img);  % 将像素值归一化到 [0, 1]
%                 img_adj = imadjust(img_normalized, [low_RGB_triplet; high_RGB_triplet], []);
%                 imwrite(img_adj, saveName);
                
                fprintf('保存至 %s \n', saveName);
    
    %             %图像对展示
    %             % figure
    %             h = montage(...
    %                 {mat2gray(FusionImg(:,:,4:-1:2)), ...
    %                 mat2gray(FusionImg(:,:,4:-1:2))}, ...
    %                 'BorderSize',10,'BackgroundColor','white')
    %             title('全色图像 (左)和多光谱图像 (右)');
    %             %}
            end
    %% 数据保存
    % 将一组数据以mat格式传到对应的文件夹里面，并保留对应的缩略图
    % 
    % saveDir =[ 'H:\Benchmark\',pathList{i}(PathLen+1:end-12)];%设置对应保存路径            
    % end
    % NumPatchs = NumPatchs+1;
    % formatSpec = '保存第%d- %d个图像对！\n';
    % fprintf(formatSpec, i, NumPatchs);
    
            %保存融合数据集重新命名
    %         Pan= patch_Pan;
    %         MS = patch_MS;
    %         MS_Up =  patch_MS_Up;
    %         Pan_LR = patch_Pan_LR;
    %         MS_LR = patch_MS_LR;
    %         MS_LR_Up = patch_MS_LR_Up;
%             output = FusionImg;
            %参数保存
%             Paras.ratio = Scale;%分辨率
%             Paras.sensor = SensorName;%传感器类型
%             Paras.intre = 'bicubic';%插值方式
         
        end
    end
        
    fprintf("已完成，请到%s查看！\n ",saveDir);
end
