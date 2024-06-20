%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
% 从FusionImg2EvaluationMetric.m中拆出来降分辨率的
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [RB,... (1) Relative Bias 相对偏差
    RV,... (2) Relative Variance 相对方差
    RSD,... (3) Relative Standard Deviation 相对标准偏差
    ...SAM,... (4) Spectral angle mapper 光谱角映射器
    RMSE_,... (5) Root Mean Square Error 均方根误差
    ERGAS_,... (6) Erreur Relative Globale Adimensionnelle de synthese 相对全局无量纲综合误差
    QAVE_,... (7) Universal Image Quality Index 通用图像质量指数
    CCMean... (8) Correlation Coefcient 相关系数
    ...mssim,... (9) Structural Similarity Index Measure 结构相似性指数测量
    ...ssim_map,... (9)
    ...SD,... (10) Standard Deviation 标准差
    ...entropy_,... (11) Entropy and cross entropy 熵和交叉熵
    ...CEMean,... (11) Entropy and cross entropy 熵和交叉熵
    ...MI,... (12) Mutual Information 互信息
    ...SFMean... (13) Spatial Frequency 空间频率
    ...QNR,... (14) (QNR)无价值参考质量指数
    ...HVS... (15) HVS Consistent Fusion Quality Assessment Index HVS一致性融合质量评估指标
    ] = FusionImg2EvaluationMetricDR(I_F,I_GT)
    

% function [D_lambda,D_S,QNRI,SAM,SCC,Q_index,SAM_index,ERGAS_index,sCC,Q2n_index,RB,RV,RSD,RMSE_,ERGAS_,QAVE_,CCMean,SD,entropy_,CEMean,MI,SFMean] = FusionImg2EvaluationMetric(I_F,I_MS,I_PAN,I_MS_Up)


%% Reference Image‑Based Quality Metrics 

%% (1) Relative Bias 相对偏差
% 相对偏差是使用输入多光谱图像的平均值与全色锐化后图像的平均值的差，再除以输入多光谱图像的平均值，理想值是0.

mean_MS = mean(I_GT(:));% 计算多光谱图像和全色图像的平均值
mean_PAN = mean(I_F(:));
RB = abs((mean_MS - mean_PAN)) / mean_MS;% 计算相对偏差

%% (2) Relative Variance 相对方差
% 输入多光谱图像和全色锐化后图像的平方差，然后将其除以参考图像平方。相对标准偏差的理想值为 0。

squared_difference = (I_GT .^2) - (I_F .^2); % 计算平方差
squared_difference = abs(squared_difference); % 计算平方差的绝对值
squared_reference = I_GT .^2; % 计算参考图像的平方

% RV = mean(squared_difference(:)) / mean(squared_reference(:)) ; % 计算相对标准偏差
RV = squared_difference ./ squared_reference ;
RV = mean(RV(:));
%% (3) Relative Standard Deviation 相对标准偏差
% 参考图像和融合多光谱波段的差值，除以参考图像的平均值，相对标准偏差的理想值为0。

RSD = abs( I_GT - I_F ) ./ mean(I_GT);
RSD = mean(RSD(:));

%% (4) Spectral angle mapper 光谱角映射器
% 多光谱图像和融合图像的矢量表示中的余弦角来表示
% .\Evluation_Metrics\Metric_Code\SAM.m
% MS = I_MS_Up; F = I_F;
% SAM = SAM(MS,F);

%% (5) Root Mean Square Error 均方根误差
% 找到参考图像和融合图像的减法，然后找到其均方根值,价值应该很低。
% .\Evluation_Metrics\Metric_Code\RMSE.m
MS = I_GT; F = I_F;
RMSE_ = RMSE(MS,F);

%% (6) Erreur Relative Globale Adimensionnelle de synthese 相对全局无量纲综合误差
% 通过RMSE 值的总和来计算。其中 p 是空间，q 是多光谱图像的光谱分辨率。 휇 i 代表平均强度，K 是多光谱带总数。
% .\Evluation_Metrics\Metric_Code\ERGAS.m
MS = I_GT; F = I_F;
ERGAS_ = ERGAS(MS,F);

%% (7) Universal Image Quality Index 通用图像质量指数qave
% 它也称为 Q 指数。其中，图像失真可以用对比度失真、亮度失真和相关性损失三个因素来表示。通用图像质量指数的取值范围在-1到+1之间。它的值必须接近1。
% 然后，对所有局部 Q 指数进行平均，以找到全局 Q 指数。 Q4 指数是此质量指标的扩展，用于四波段数据集。如果有四个以上的光谱带，则 Q4 的广义形式是 Q2n 指数。
% .\Evluation_Metrics\Metric_Code\QAVE.m
MS = I_GT; PANMS = I_F;
QAVE_ = QAVE(MS,PANMS);

%% (8) Correlation Coefcient 相关系数
% 该指标可找出全色锐化图像和参考图像之间的相关性。相关系数计算为：
% 其中 fm 和 rm 分别是融合多光谱波段和参考多光谱波段，大小为 M × N。此外，− rm 和 − fm 是平均值。
% .\Evluation_Metrics\spectral_metric\correlation_coefficient.m

% 循环处理每个波段
[num_bands] = size(I_F,3);% 获取数组的维度信息
for band_idx = 1:num_bands
    % 从第一个数组中提取单波段图像
    img1 = I_GT(:, :, band_idx);

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
% img1 = I_MS; img2 = I_F;
% [mssim, ssim_map] = SSIM_index(img1, img2, K, window, L);





end