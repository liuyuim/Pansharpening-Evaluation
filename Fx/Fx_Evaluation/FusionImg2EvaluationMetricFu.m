%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
% 从FusionImg2EvaluationMetric.m中拆出来全分辨率的
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [...RB,... (1) Relative Bias 相对偏差
    ...RV,... (2) Relative Variance 相对方差
    ...RSD,... (3) Relative Standard Deviation 相对标准偏差
    ...SAM,... (4) Spectral angle mapper 光谱角映射器
    ...RMSE_,... (5) Root Mean Square Error 均方根误差
    ...ERGAS_,... (6) Erreur Relative Globale Adimensionnelle de synthese 相对全局无量纲综合误差
    ...QAVE_,... (7) Universal Image Quality Index 通用图像质量指数
    ...CCMean,... (8) Correlation Coefcient 相关系数
    ...mssim,... (9) Structural Similarity Index Measure 结构相似性指数测量
    ...ssim_map,... (9)
    SD,... (10) Standard Deviation 标准差
    entropy_,... (11) Entropy and cross entropy 熵和交叉熵
    CEMean,... (11) Entropy and cross entropy 熵和交叉熵
    ...MI,... (12) Mutual Information 互信息
    SFMean... (13) Spatial Frequency 空间频率
    ...QNR,... (14) (QNR)无价值参考质量指数
    ...HVS... (15) HVS Consistent Fusion Quality Assessment Index HVS一致性融合质量评估指标
    ] = FusionImg2EvaluationMetricFu(I_F,I_MS_Up)
    



%% 5.1 Quality metrics without reference image
%% (10) Standard Deviation 标准差
% 它测量融合图像中的对比度，其中对比度值越高表示信息丰富度。该措施对于无噪声图像的情况是有效的。
% .\Evluation_Metrics\spatial_metric\Standard_deviation.m
img = I_F;
SD = Standard_deviation(img);

%% (11) Entropy and cross entropy 熵和交叉熵
% 熵和交叉熵它再次衡量图像的丰富度。熵的计算公式为：
% 其中要分析的图像的动态范围由 K 表示。第k个灰度级出现概率用p(t)表示。交叉熵发现原始图像和融合图像之间的相似性。这两项措施都应该具有较高的价值，以获得更好的结果。
% .\Evluation_Metrics\spectral_metric\cross_entropy.m

% Entropy熵
entropy_ = entropy_1(I_F);

% cross entropy交叉熵
% 循环处理每个波段
[num_bands] = size(I_F,3);% 获取数组的维度信息
for band_idx = 1:num_bands
    % 从第一个数组中提取单波段图像 %cross entropy交叉熵暂时借用一下I_MS_Up
    img1 = I_MS_Up(:, :, band_idx);

    % 从第二个数组中提取单波段图像
    img2 = I_F(:, :, band_idx);

    % 在这里执行您想要针对单波段图像执行的操作，例如保存、显示等
    CE(1,band_idx) = cross_entropy(img1,img2);
end
CE1=CE(1,1); CE2=CE(1,2); CE3=CE(1,3); CE4=CE(1,4); CEMean=(CE1+CE2+CE3+CE4)/4;


%% (12) Mutual Information 互信息
% 通过互信息计算联合概率分布。它是在参考图像和融合图像之间计算的。计算互信息的方程为：
% 这里，I 和 F 是随机变量，分别代表源图像和融合图像。 PI(i) 和 PF(f) 是概率分布。

% % 计算图像的灰度级数（假设是8位图像）
% num_gray_levels = 256;
% % 计算参考图像和融合图像的联合直方图
% joint_histogram = histcounts2(I(:), F(:), num_gray_levels, num_gray_levels);
% % 计算联合概率分布
% joint_prob = joint_histogram / sum(joint_histogram(:));
% % 计算参考图像的概率分布
% prob_I = sum(joint_prob, 2);
% % 计算融合图像的概率分布
% prob_F = sum(joint_prob, 1);
% % 计算互信息
% mutual_information = 0;
% for i = 1:num_gray_levels
%     for f = 1:num_gray_levels
%         if joint_prob(i, f) > 0
%             mutual_information = mutual_information + joint_prob(i, f) * log2(joint_prob(i, f) / (prob_I(i) * prob_F(f)));
%         end
%     end
% end


%% (13) Spatial Frequency 空间频率
% 通过空间频率计算总体活动水平,它可以用列频率和行频率来表示。空间测量必须具有较高的价值才能获得更好的结果,越大越好。，RF 是行频率。 CF是列频率。
% .\Evluation_Metrics\spatial_metric\spatial_frequency.m

% 循环处理每个波段
[num_bands] = size(I_F,3);% 获取数组的维度信息
for band_idx = 1:num_bands
    % 从第一个数组中提取单波段图像    
    img = I_F(:, :, band_idx);

    % 在这里执行您想要针对单波段图像执行的操作，例如保存、显示等
    SF(1,band_idx) = spatial_frequency(img);
end
SF1=SF(1,1); SF2=SF(1,2); SF3=SF(1,3); SF4=SF(1,4); SFMean=(SF1+SF2+SF3+SF4)/4;

%% (14) Quality with no value reference (QNR) Index 无价值参考质量(QNR)指数
% 包括空间和光谱失真指数，A and B are weight coefcients.
% .\Evluation_Metrics\spectral_metric\QNR.m 这个是空的

%% (15) HVS Consistent Fusion Quality Assessment Index HVS一致性融合质量评估指标



end