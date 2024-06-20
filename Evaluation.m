%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 【5 融合结果的批量评价】
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 目视观察
% >> size(ms)ans =   100    16    16     8
% K>> imshow(mat2gray(I_MS_LR(:,:,4:-1:2)))
% K>> imwrite(mat2gray(I_MS_LR(:,:,4:-1:2)),'1.bmp')   视频里的两条命令
% imshow(I_MS_LR(1:256,1:256,1),[])
% 
% >> imwrite(mat2gray(MS(:,:,4:-1:2)),'1_MS.bmp')
% >> imwrite(mat2gray(Pan(:,:,1)),'1_Pan.bmp')

%% Bmp 颜色有点问题，能用来批量目视用
clc;clear;close all;addpath(genpath('.\Fx\'));

% 指定目录
% FusionImgYijiPath ='F:\Demo\Data_Evaluation\AmodelOutput_Fu128\5AmodelOutput_Fu128'; %测试集经过深度学习test代码得出的五个假设融合结果，
% saveDir = 'F:\Demo\Data_Evaluation\AmodelOutput_Fu128\5AmodelOutput_Fu128_bmp';  %设置对应保存路径 Fu评价用Fu数据集
% MatOutputErjiDir2Bmp (FusionImgYijiPath,Step,saveDir) %mat里面是output

%% 指定目录
clc;clear;close all;addpath(genpath('.\Fx\'));
Step = 1;
FusionImgYijiPath =''; %测试集经过深度学习test代码得出的五个假设融合结果，
saveDir = '';  %设置对应保存路径 Fu评价用Fu数据集
MatMSErjiDir2Bmp (FusionImgYijiPath,Step,saveDir)

%% 循环
clc;clear;close all;addpath(genpath('.\Fx\'));
Step = 5; NetNames = {'WSDFNet'};  % 'PanNet','LPPN','WSDFNet'
for i = 1:numel(NetNames)
    NetName = NetNames{i};
    SensorNames = {'GF1','QB'}; %{'GF1','IK','QB','WV2','WV3','WV4'} {'GF1','GF2','JL1','QB','WV2','WV3'}  'GF1','GF2','IK','JL1','QB','WV2','WV3','WV4'
    for j = 1:numel(SensorNames)
    Sensor_Data = strcat(SensorNames{j}, '_Data');    % 或者    Sensor_Data = SensorNames{i} + "_Data";
    Sensor_Net = strcat(SensorNames{j},'_',NetName); 
        for Sizes = [64,32] % 1024,512,256,128,64,32
            Size = num2str(Sizes);  % AmodelOutput_Fu
                        
            HypothesisOutput_Fu = ['HypothesisOutput_Fu',Size]; HypothesisOutput_Fu_bmp = ['HypothesisOutput_Fu',Size,'_bmpHsi']; 
            
            FusionImgYijiPath = fullfile('F:\Demo\Data_Evaluation',Sensor_Net,HypothesisOutput_Fu); %测试集经过深度学习test代码得出的五个假设融合结果，
            
            saveDir = fullfile('F:\Demo\Data_Evaluation',Sensor_Net,HypothesisOutput_Fu_bmp);   
            
            MatOutputErjiDir2Bmp (FusionImgYijiPath,Step,saveDir)
        end
    end
end

%% Tif (用于后续ENVI中处理)
clc;clear;close all;addpath(genpath('.\Fx\'));

Sensor = 'WV3';    Net = 'WSDFNet'; 
InputMatFileName = ['F:\Demo\Data_Evaluation\',Sensor,'_',Net,'\HypothesisOutput_Fu1024\HypothesisIn','WV3','model\27.mat'];
% InputMatFileName = ['F:\Demo\Data_Evaluation\',Sensor,'_',Net,'\HypothesisOutput_Fu1024\HypothesisIn',Sensor,'model\20.mat'];
OutputTifFilename = ['F:\Demo\Data_Evaluation\',Sensor,'_',Net,'\27.tif'];
Mat2Tif(InputMatFileName,OutputTifFilename);

%% 量化评价
% 采用DeepLearning数据集进行融合实验，在高低两个分辨率尺度上对融合结果进行评价
% 对比方法采用 Pansharpening Tool ver 1.3中的方法；
% 对比指标是 Fu：D_lambda,D_S,QNRI,SAM,SCC。DR：Q_index, SAM_index, ERGAS_index, sCC, Q2n_index。
% Full-resolution：
% (QNR)quality with no reference 无参考质量 The bigger the better  "1"
% Dm in the QNR assesses the spectral distortion, 光谱失真，The smaller the better  "0"
% Ds inthe QNR assesses the spatial distortion ,空间失真，The smaller the better  "0"
% sCC spatial Correlation Coefficient . 融合图像和 PAN 图像之间的空间相关系数 "1"
% Reduced-resolution：
% (ERGAS) Erreur Relative Globale Adimensionnelle de Synthèse 无量纲全局相对综合误差  "0"
% (SAM) spectra mapper angle 光谱映射器角度 计算融合图像和参考图像中相应像素之间的角度  "0"
% Q 指数通用图像质量指数 "1"
% Q2n index Q4 指数的推广适用于评估具有大于四个光谱带的图像 The bigger the better
% (RMSE)均方根误差 "0"

% 对多种传感器评价后生成的"5种评价指标*100张图片的"二维数值统计表MatrixResults，再次罗列起来生成MatrixAll.mat 
% MatrixAll.mat 
% 第一维是若干种假设融合图像的评价结果，HypothesInGF1,HypothesInGF2...
% 第二维是五个评价指标,如 D_lambda,D_S,QNRI,SAM,SCC
% 第三维是100张图片
% 
% format short g 为默认格式，保留四位小数，是四舍五入；
% format long g  长精度，保留16位小数，是四舍五入；
% format bank g  保留两位小数，不是四舍五入；
% format rat g   采用分数或者整数的形式表示结果。

clc;clear;close all;addpath(genpath('.\Fx\'));


% Full-resolution
% FusionImgYijiPath='F:\Demo\Data_Evaluation\GF1_WSDFNet\HypothesisOutput_Fu1024'; %里面得套个二级目录 测试集经过深度学习test代码得出的五个假设融合结果，
% BenchYijiPath='F:\Demo\Data_Evaluation\GF1_Data\Test_Fu1024'; %测试集
% saveDir = 'F:\Demo\Data_Evaluation\GF1_WSDFNet\EvaluateFu1024';  %设置对应保存路径 Fu评价用Fu数据集
% FusionImg2EvaluationFu (FusionImgYijiPath,BenchYijiPath,saveDir);

NetNames = {'WSDFNet'}; %'PanNet','LPPN','WSDFNet'
for i = 1:numel(NetNames)
    NetName = NetNames{i};
    SensorNames = {'GF1','QB'}; %{'GF1','IK','QB','WV2','WV3','WV4'} {'GF1','GF2','JL1','QB','WV2','WV3'}  'GF1','GF2','IK','JL1','QB','WV2','WV3','WV4'
    for j = 1:numel(SensorNames)
    Sensor_Data = strcat(SensorNames{j}, '_Data');    % 或者    Sensor_Data = SensorNames{i} + "_Data";
    Sensor_Net = strcat(SensorNames{j},'_',NetName); 
        for Sizes = [64,32] % 1024,512,256,128,64,32
            Size = num2str(Sizes);  
            HypothesisOutput_Fu = ['HypothesisOutput_Fu',Size]; Test_Fu = ['Test_Fu',Size]; Evaluate_Fu = ['EvaluateFu',Size];

            FusionImgYijiPath = fullfile('F:\Demo\Data_Evaluation',Sensor_Net,HypothesisOutput_Fu); %测试集经过深度学习test代码得出的五个假设融合结果，
            TestYijiPath = fullfile('F:\Demo\Data_Evaluation',Sensor_Data,Test_Fu); %测试集
            saveDir = fullfile('F:\Demo\Data_Evaluation',Sensor_Net,Evaluate_Fu);  %设置对应保存路径 注意，Fu评价要用Fu数据集
            FusionImg2EvaluationFu (FusionImgYijiPath,TestYijiPath,saveDir);

        end
    end
end


% Reduced-resolution 
% FusionImgYijiPath='G:\AFusionGroup\Shiyan\shiyan20230903\GF1_WSDFNet\HypothesisOutput_DR'; %测试集经过深度学习test代码得出的五个假设融合结果，
% BenchYijiPath='G:\AFusionGroup\Shiyan\shiyan20230903\GF1_Data\Benchmark\'; %Benchmark作为测试集
% saveDir = 'G:\AFusionGroup\Shiyan\shiyan20230903\GF1_WSDFNet\EvaluateDR';  %设置对应保存路径 
% FusionImg2EvaluationDR (FusionImgYijiPath,BenchYijiPath,saveDir);
NetNames = {'WSDFNet'}; %'PanNet','LPPN','WSDFNet'
for i = 1:numel(NetNames)
    NetName = NetNames{i};
    SensorNames = {'GF1','QB'}; %{'GF1','IK','QB','WV2','WV3','WV4'} {'GF1','GF2','JL1','QB','WV2','WV3'}  'GF1','GF2','IK','JL1','QB','WV2','WV3','WV4'
    for j = 1:numel(SensorNames)
    Sensor_Data = strcat(SensorNames{j}, '_Data');    % 或者    Sensor_Data = SensorNames{i} + "_Data";
    Sensor_Net = strcat(SensorNames{j},'_',NetName); 
        for Sizes = [64,32] % 1024,512,256,128,64,32
            Size = num2str(Sizes);                        
            % HypothesisOutput_Fu = ['AmodelOutput_Fu',Size]; Test_Fu = ['Test_Fu',Size]; Evaluate_Fu = ['AmodelOutputEvaluate_Fu',Size];
            
            HypothesisOutput_Fu = ['HypothesisOutput_DR',Size]; Test_Fu = ['Test_DR',Size]; Evaluate_Fu = ['EvaluateDR',Size];

            FusionImgYijiPath = fullfile('F:\Demo\Data_Evaluation',Sensor_Net,HypothesisOutput_Fu); %测试集经过深度学习test代码得出的五个假设融合结果，
            TestYijiPath = fullfile('F:\Demo\Data_Evaluation',Sensor_Data,Test_Fu); %测试集
            saveDir = fullfile('F:\Demo\Data_Evaluation',Sensor_Net,Evaluate_Fu);  %设置对应保存路径 注意，Fu评价要用Fu数据集
            FusionImg2EvaluationDR (FusionImgYijiPath,TestYijiPath,saveDir);
            
        end
    end
end
fprintf("所有 融合影像评价 完成，脚本程序结束！\n");


%% 5.3 量化评价大全 全分辨率+降分辨率
% 对比方法采用 Pansharpening Tool ver 1.3中的函数 + Evluation_Metrics 中的函数 
clc;clear;close all;addpath(genpath('.\Fx\'));

NetNames = {'WSDFNet'}; %'PanNet','LPPN','WSDFNet'
for i = 1:numel(NetNames)
    NetName = NetNames{i};
    SensorNames = {'GF1','GF2','JL1','QB','WV2','WV3'}; %{'GF1','IK','QB','WV2','WV3','WV4'} {'GF1','GF2','JL1','QB','WV2','WV3'}  'GF1','GF2','IK','JL1','QB','WV2','WV3','WV4'
    for j = 1:numel(SensorNames)
    Sensor_Data = strcat(SensorNames{j}, '_Data');    % 或者    Sensor_Data = SensorNames{i} + "_Data";
    Sensor_Net = strcat(SensorNames{j},'_',NetName); 
        for Sizes = [256,128,64,32] % 1024,512,256,128,64,32
            Size = num2str(Sizes);                        
%             HypothesisOutput_Fu = ['AmodelOutput_Fu',Size]; Test_Fu = ['Test_Fu',Size]; Evaluate_Fu = ['AmodelOutputEvaluate_Fu',Size];
            % HypothesisOutput_Fu = ['HypothesisOutput_Fu',Size]; Test_Fu = ['Test_Fu',Size]; Evaluate_Fu = ['Evaluate_Fu',Size];
            HypothesisOutput_Fu = ['HypothesisOutput_DR',Size]; Test_Fu = ['Test_DR',Size]; Evaluate_Fu = ['Evaluate_DR',Size];

            FusionImgYijiPath = fullfile('F:\Demo\Data_Evaluation',Sensor_Net,HypothesisOutput_Fu); %测试集经过深度学习test代码得出的五个假设融合结果，
            TestYijiPath = fullfile('F:\Demo\Data_Evaluation',Sensor_Data,Test_Fu); %测试集
            saveDir = fullfile('F:\Demo\Data_Evaluation',Sensor_Net,Evaluate_Fu);  %设置对应保存路径 注意，Fu评价要用Fu数据集
            FusionImg2EvaluationAll (FusionImgYijiPath,TestYijiPath,saveDir);
            
        end
    end
end

fprintf("所有 融合影像评价 完成，脚本程序结束！\n");

%% 5.4 评价结果自助打包
clc; clear; close all; addpath(genpath('.\Fx\'));

NumImgStart = 1; NumImgEnd = 100; % 文件夹中按自然数排序第多少张的范围，不是文件名范围
for Sizes = [64]  %% Size,512,256,128,64,32
    Size = num2str(Sizes);
    
    NetNames = {'WSDFNet'}; %'PanNet','LPPN','WSDFNet'
    for i = 1:numel(NetNames)
        NetName = NetNames{i};
        
        SensorNames = {'GF1','QB'}; %% Sensor 'GF1','GF2','IK','JL1','QB','WV2','WV3','WV4'
        for j = 1:numel(SensorNames)
            Sensor_Net = strcat(SensorNames{j},'_',NetName); % 或者    Sensor_Data = SensorNames{i} + "_Data";Sensor_Data = strcat(SensorNames{i}, '_Data');  
            Evaluate_Fu = ['EvaluateFu',Size];
            EvaluationDir = fullfile('F:\Demo\Data_Evaluation',Sensor_Net,Evaluate_Fu);
            saveDirName = [Evaluate_Fu,'_GF1-QB']; % ['EvaluateDRFu256_GF1-QB-WV4_',Size];
            Evaluation2RepackMatrix (SensorNames,EvaluationDir,saveDirName,NumImgStart,NumImgEnd) %'F:\Demo\Data_Evaluation\IK_WSDFNet\Evaluate_Fu128',...
          
        end
    end
end

