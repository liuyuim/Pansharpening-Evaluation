%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           Full resolution quality indexes. 
% 
% Interface:
%           [D_lambda,D_S,QNR_index,SAM_index,sCC] = indexes_evaluation_FS(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,sensor,tag,ratio)
%
% Inputs:
%           I_F:                Fused Image;
%           I_MS_LR:            MS image;
%           I_PAN:              Panchromatic image;
%           L:                  Image radiometric resolution; 
%           th_values:          Flag. If th_values == 1, apply an hard threshold to the dynamic range;
%           I_MS:               MS image upsampled to the PAN size;
%           sensor:             String for type of sensor (e.g. 'WV2','IKONOS');
%           tag:                Image tag. Often equal to the field sensor. It makes sense when sensor is 'none'. It indicates the band number;
%           ratio:              Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           D_lambda:           D_lambda index;
%           D_S:                D_S index;
%           QNR_index:          QNR index;
%           SAM_index:          Spectral Angle Mapper (SAM) index between fused and MS image;
%           sCC:                spatial Correlation Coefficient between fused and PAN images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [RB,... (1) Relative Bias 相对偏差
    RV,... (2) Relative Variance 相对方差
    RSD,... (3) Relative Standard Deviation 相对标准偏差
    ...SAM,... (4) Spectral angle mapper 光谱角映射器
    RMSE_,... (5) Root Mean Square Error 均方根误差
    ERGAS_,... (6) Erreur Relative Globale Adimensionnelle de synthese 相对全局无量纲综合误差
    QAVE_,... (7) Universal Image Quality Index 通用图像质量指数
    CCMean,... (8) Correlation Coefcient 相关系数
    ...mssim,... (9) Structural Similarity Index Measure 结构相似性指数测量
    ...ssim_map,... (9)
    SD,... (10) Standard Deviation 标准差
    entropy_,... (11) Entropy and cross entropy 熵和交叉熵
    CEMean,... (11) Entropy and cross entropy 熵和交叉熵
    ...MI,... (12) Mutual Information 互信息
    SFMean... (13) Spatial Frequency 空间频率
    ...QNR,... (14) (QNR)无价值参考质量指数
    ...HVS... (15) HVS Consistent Fusion Quality Assessment Index HVS一致性融合质量评估指标
    ] = FusionImg2EvaluationMetric(I_F,I_MS,I_PAN,I_MS_Up)
    

% function [D_lambda,D_S,QNRI,SAM,SCC,Q_index,SAM_index,ERGAS_index,sCC,Q2n_index,RB,RV,RSD,RMSE_,ERGAS_,QAVE_,CCMean,SD,entropy_,CEMean,MI,SFMean] = FusionImg2EvaluationMetric(I_F,I_MS,I_PAN,I_MS_Up)


%% Reference Image‑Based Quality Metrics % 此处I_MS_Up当gt

%% (1) Relative Bias 相对偏差
% 相对偏差是使用输入多光谱图像的平均值与全色锐化后图像的平均值的差，再除以输入多光谱图像的平均值，理想值是0.

mean_MS = mean(I_MS_Up(:));% 计算多光谱图像和全色图像的平均值
mean_PAN = mean(I_F(:));
RB = abs((mean_MS - mean_PAN)) / mean_MS;% 计算相对偏差

%% (2) Relative Variance 相对方差
% 输入多光谱图像和全色锐化后图像的平方差，然后将其除以参考图像平方。相对标准偏差的理想值为 0。

squared_difference = (I_MS_Up .^2) - (I_F .^2); % 计算平方差
squared_difference = abs(squared_difference); % 计算平方差的绝对值
squared_reference = I_MS_Up .^2; % 计算参考图像的平方

% RV = mean(squared_difference(:)) / mean(squared_reference(:)) ; % 计算相对标准偏差
RV = squared_difference ./ squared_reference ;
RV = mean(RV(:));
%% (3) Relative Standard Deviation 相对标准偏差
% 参考图像和融合多光谱波段的差值，除以参考图像的平均值，相对标准偏差的理想值为0。

RSD = abs( I_MS_Up - I_F ) ./ mean(I_MS_Up);
RSD = mean(RSD(:));

%% (4) Spectral angle mapper 光谱角映射器
% 多光谱图像和融合图像的矢量表示中的余弦角来表示
% .\Evluation_Metrics\Metric_Code\SAM.m
% MS = I_MS_Up; F = I_F;
% SAM = SAM(MS,F);

%% (5) Root Mean Square Error 均方根误差
% 找到参考图像和融合图像的减法，然后找到其均方根值,价值应该很低。
% .\Evluation_Metrics\Metric_Code\RMSE.m
MS = I_MS_Up; F = I_F;
RMSE_ = RMSE(MS,F);

%% (6) Erreur Relative Globale Adimensionnelle de synthese 相对全局无量纲综合误差
% 通过RMSE 值的总和来计算。其中 p 是空间，q 是多光谱图像的光谱分辨率。 휇 i 代表平均强度，K 是多光谱带总数。
% .\Evluation_Metrics\Metric_Code\ERGAS.m
MS = I_MS_Up; F = I_F;
ERGAS_ = ERGAS(MS,F);

%% (7) Universal Image Quality Index 通用图像质量指数qave
% 它也称为 Q 指数。其中，图像失真可以用对比度失真、亮度失真和相关性损失三个因素来表示。通用图像质量指数的取值范围在-1到+1之间。它的值必须接近1。
% 然后，对所有局部 Q 指数进行平均，以找到全局 Q 指数。 Q4 指数是此质量指标的扩展，用于四波段数据集。如果有四个以上的光谱带，则 Q4 的广义形式是 Q2n 指数。
% .\Evluation_Metrics\Metric_Code\QAVE.m
MS = I_MS_Up; PANMS = I_F;
QAVE_ = QAVE(MS,PANMS);

%% (8) Correlation Coefcient 相关系数
% 该指标可找出全色锐化图像和参考图像之间的相关性。相关系数计算为：
% 其中 fm 和 rm 分别是融合多光谱波段和参考多光谱波段，大小为 M × N。此外，− rm 和 − fm 是平均值。
% .\Evluation_Metrics\spectral_metric\correlation_coefficient.m

% 循环处理每个波段
[num_bands] = size(I_F,3);% 获取数组的维度信息
for band_idx = 1:num_bands
    % 从第一个数组中提取单波段图像
    img1 = I_MS_Up(:, :, band_idx);

    % 从第二个数组中提取单波段图像
    img2 = I_F(:, :, band_idx);

    % 在这里执行您想要针对单波段图像执行的操作，例如保存、显示等
    CC(1,band_idx) = correlation_coefficient(img1,img2);
end
CC1=CC(1,1); CC2=CC(1,2); CC3=CC(1,3); CC4=CC(1,4); CCMean=(CC1+CC2+CC3+CC4)/4;

% 方案2：将图像reshape成矩阵形式，每行表示一个像素，每列表示一个波段
% 
% image1_matrix = reshape(I_MS_Up, [], 1);
% image2_matrix = reshape(I_F, [], 1);
% % 计算相关系数矩阵
% correlation_matrix = corr(image1_matrix, image2_matrix);
% CC = correlation_coefficient(image1_matrix,image2_matrix);
% % 显示相关系数矩阵
% imagesc(correlation_matrix); colorbar; title('Correlation Coefficient Matrix');

%% (9) Structural Similarity Index Measure 结构相似性指数测量
% 融合图像和参考图像可能具有通过此测量计算的结构相似性。人类可以直观地看到原始图像和融合图像中结构细节的差异。这种差异说明了融合性能。该度量的范围为 [-1,1]，其中 1 表示图像相同。
% .\Evluation_Metrics\spatial_metric\SSIM_index.m
K = [0.01 0.03]; % K = [0.05 0.05];
window = fspecial('gaussian', 11, 1.5); % window = ones(8);
L = 255; % L = 100;
img1 = I_MS; img2 = I_F;
% [mssim, ssim_map] = SSIM_index(img1, img2, K, window, L);

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
    % 从第一个数组中提取单波段图像
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