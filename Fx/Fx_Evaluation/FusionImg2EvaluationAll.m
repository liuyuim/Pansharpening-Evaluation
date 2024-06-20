%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = FusionImg2EvaluationAll (FusionImgYijiPath,TestYijiPath,saveDir)

    % TestData是不变的，利用它得到NumImgs
    % BenchErjiPath = fullfile(TestYijiPath); % Bench 只有一级
    Test_list = dir([TestYijiPath,'\','*.mat']) ;
    NumImgs = size(Test_list,1);

    % 用FusionImgYijiPath 遍历出二级目录名，共用 GF1_GengDi GF1_LinDi   GF1_WeiBiaoDuoLei ...
    ErjiDir_list = dir(FusionImgYijiPath) ;  % 二级目录列表
    ErjiDir_list_Nums = size(ErjiDir_list,1);  % 二级目录个数 包括 .和..
    
    % 定义最后保存所有东西的矩阵
    Matrix_Fu = []; 

    for i_ErjiDir = 3 : ErjiDir_list_Nums
        %列出当前二级文件夹内所有的mat
        FusionImgErjiPath = fullfile(FusionImgYijiPath,ErjiDir_list(i_ErjiDir).name); %.\DataDL_PannetOutput\GF1_GengDi
        FusionImg_list = dir([FusionImgErjiPath,'\','*.mat']) ;
        % 设置xlsx文件名
        XlsxName = strcat("EvaluateReport", string(datetime, 'yyyy-MM-dd-HH-mm-ss'), '.xlsx');
        
        % 在当前二级目录处理每一个mat
        for i_NumImgs = 1:NumImgs
        
            formatSpec = '开始处理二级目录 %s！%d个图像中第%d个！... \n';
            fprintf(formatSpec,ErjiDir_list(i_ErjiDir).name, NumImgs, i_NumImgs);
    
            % 校验 当前从 Output和TestData文件夹 分别取出的 mat文件名 是否一致
            %验证两者是否一致
            if ~isequal(FusionImg_list(i_NumImgs).name, Test_list(i_NumImgs).name)
                fprintf("当前从 Output和TestData文件夹分别取出的 mat文件名 不一致");
                break;
            end
            
            % 然后再正常运行
            
            %把mat文件加载进来
            FusionImgPath = [FusionImg_list(i_NumImgs).folder,'\',FusionImg_list(i_NumImgs).name]; %TestOutput_list列表中的第i个目录和文件名拼成要加载的mat路径 如 E:\LiuYu\FusionEvaluateExperiment\DataDL_PannetOutput\j1p1.mat
            FusionImgDate = load(FusionImgPath); 
            TestDataPath = [Test_list(i_NumImgs).folder,'\',Test_list(i_NumImgs).name]; %TestData_list列表中的第i个目录和文件名拼成要加载的mat路径  E:\LiuYu\FusionEvaluateExperiment\DataDLPre_1TestData\j1p1.mat
            TestData = load(TestDataPath); 
            
                         
            
        
            
          
            %% 逆归一化
                        
            % 归一化，逆归一化原理：
            % 1、默认的归一化范围是（-1，1），使用mapminmax(data,0,1)将范围控制在（0，1）。
            % 2、按行归一化，矩阵则每行归一化一次。若要完全归一化，则
            % I_F = double(FusionImgDate.output);
            % FlattenedData = I_F(:)'; % 展开矩阵为一列，然后转置为一行。
            % MappedFlattened = mapminmax(FlattenedData, 0, 1); % 归一化。
            % MappedFlattened = MappedFlattened*2047 ;
            % I_F = reshape(MappedFlattened, size(I_F)); %还原为原始矩阵形式。此处不需转置回去，因为reshape恰好是按列重新排序
            % 
            % 逆归一化方案1： I_F = (mat2gray(FusionImgDate.output))*(2^11-1); %I_F: Fused Image;  mat2gray后即转为double
            % 逆归一化方案2： I_F = (double(FusionImgDate.output))*(2^11-1); 发现这个效果好一点用的这个
            
            if contains(ErjiDir_list(i_ErjiDir).name,"GF1"|"GF2") %contains确定字符串中是否有模式,matches确定模式是否与字符串匹配
                I_F = (double(FusionImgDate.output))*(2^10-1);
                fprintf("当前HypothesisDir 匹配到是GF1/GF2,使用*1023！...");
            elseif contains(ErjiDir_list(i_ErjiDir).name,"A"|"IK"|"JL"|"QB"|"WV")
                I_F = (double(FusionImgDate.output))*(2^11-1);
                fprintf("当前HypothesisDir 匹配到是A/IK/JL/QB/WV,使用*2047！...");
            else
                fprintf("当前从HypothesisDir 匹配不到是哪种传感器的假设文件夹，检查代码Fx/FusionImg2Evaluation.m");
                break;
            end

            cd  '.\Toolbox\Pansharpening Tool ver 1.3\'

            % indexes_evaluation_FS和indexes_evaluation 用的a critical代码，变量命名按他的体系
            %% indexes_evaluation_FS
            I_MS_LR = double(TestData.ms); % MS image;
            I_MS =  double(TestData.lms); % MS image upsampled to the PAN size;
            I_PAN = double(TestData.pan); %Pan
            
            Params = TestData.Paras;
            
            
            % Threshold values out of dynamic range
            thvalues = 0;
            L = ceil(log2(double(max(I_PAN(:)))+1));% Radiometric Resolution
            sensor = Params.sensor;
            im_tag =  Params.sensor;
            ratio = Params.ratio;
            
            %     Params = imgData.Paras;
            %     Paras.ratio = Scale;%分辨率
            %     Paras.sensor = SensorName;%传感器类型
            %     Paras.intre = 'bicubic';%插值方式
            
            t2=tic;           

%           function [D_lambda,D_S,QNR_index,SAM_index,sCC] = indexes_evaluation_FS(I_F,I_MS_LR,I_PAN,L,th_values,I_MS,sensor,tag,ratio)
            [D_lambda,D_S,QNRI,SAM,SCC] = indexes_evaluation_FS(I_F,I_MS_LR,I_PAN,L,thvalues,I_MS,sensor,im_tag,ratio);
            
            %% indexes_evaluation
            I_GT = double(TestData.lms); %ground truth 这里indexes_evaluation只是借用lms当I_GT，具体按需求来。

            Qblocks_size = 32;
            flag_cut_bounds = 0;%不进行裁切
            dim_cut = 0;%裁剪的大小不设置

            [Q_index, SAM_index, ERGAS_index, sCC, Q2n_index] = indexes_evaluation(I_F,I_GT,ratio,L,Qblocks_size,flag_cut_bounds,dim_cut,thvalues);
            
            cd ../../ %返回主文件夹,从 '.\Toolbox\Pansharpening Tool ver 1.3\' 返回
            cd  '.\Fx\Evluation_Metrics\'

            %% FusionImg2EvaluationMetric
            % 这里变量命名按benchmark的体系
            I_MS = double(TestData.ms); % MS image;
            I_MS_Up =  double(TestData.lms); % MS image upsampled to the PAN size;
            
            [RB,... (1) Relative Bias 相对偏差
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
            ] = FusionImg2EvaluationMetric(I_F,I_MS,I_PAN,I_MS_Up);
            %% 写入表中
            MatrixResult_Fu(1,:) = [D_lambda,D_S,QNRI,SAM,SCC,...  
                Q_index,SAM_index,ERGAS_index,sCC,Q2n_index,...
                RB,RV,RSD,RMSE_,ERGAS_,QAVE_,CCMean,SD,entropy_,CEMean,SFMean];
            
            cd ../../ %返回主文件夹，从 '\Evluation_Metrics\' 返回

            time_=toc(t2);
            fprintf('Elaboration time Fu: %.2f [sec]\n',time_);
        
            
%           MatrixResults_Fu(:,:,i)= cat(1, MatrixResults_Fu, MatrixResult_Fu);;%利用cat联结（按第几维来联结，被联结的图，要联结的）
            MatrixResults_Fu(1,:,i_NumImgs) = MatrixResult_Fu;
         
            % 保存每组图像的融合结果
            
%             saveName = fullfile(saveDir,ErjiDir_list(ErjiDir_i).name,[num2str(i),'.mat']);
            
            saveErjiDir_Path = fullfile(saveDir,ErjiDir_list(i_ErjiDir).name,['\']); %.\DataDL_PannetEvaluate\GF1_GengDi
            if ~exist(saveErjiDir_Path,'dir')%待保存的图像文件夹不存在，就建文件夹
                mkdir(saveErjiDir_Path)            
            end

            saveName = fullfile(saveDir,ErjiDir_list(i_ErjiDir).name,FusionImg_list(i_NumImgs).name);
            save(saveName, 'FusionImgPath','MatrixResult_Fu');
            
            %输出到xlsx https://ww2.mathworks.cn/help/matlab/ref/writematrix.html
            
            saveXlsxName = fullfile(saveDir,ErjiDir_list(i_ErjiDir).name,XlsxName);
            XlsxTitle = [ "D_lambda", "D_S", "QNRI", "SAM", "SCC", ...
                "Q_index", "SAM_index", "ERGAS_index", "sCC", "Q2n_index", ...
                "RB", "RV", "RSD", "RMSE_", "ERGAS_", "QAVE_", "CCMean", "SD", "entropy_", "CEMean", "SFMean", "Path"]; % 指标标题
%             xlswrite(saveXlsxName, XlsxTitle,['sheet',i_ErjiDir-2],'A1');
            
            writematrix(XlsxTitle,saveXlsxName)
            writematrix(MatrixResult_Fu,saveXlsxName,'WriteMode','append')
            
%             XlsxValue = cat(2,FusionImgPath,mat2str(MatrixResult_Fu));
%             XlsxValue = [FusionImgPath,num2str(MatrixResult_Fu)];
%             writematrix(XlsxValue,saveXlsxName,'WriteMode','append')

            writematrix(FusionImgPath,saveXlsxName,'Range',['x',num2str(i_NumImgs+1)]);   % F
            
            formatSpec = '已将当前图片评价结果mat文件保存至目录 %s！\n';
            fprintf(formatSpec,saveName);
        end
        
        %开始统计 
    
        %计算均值
        Mean_Fu = mean(MatrixResults_Fu,3);
        writematrix('均值',saveXlsxName,'WriteMode','append')
        writematrix(Mean_Fu,saveXlsxName,'WriteMode','append')        
        %计算中值 
        median_Fu = median(MatrixResults_Fu,3);
        writematrix('中值',saveXlsxName,'WriteMode','append')
        writematrix(median_Fu,saveXlsxName,'WriteMode','append')             
        %计算最大元素和最小元素
        max_Fu = max(MatrixResults_Fu,[],3);    
        writematrix('最大值',saveXlsxName,'WriteMode','append')
        writematrix(max_Fu,saveXlsxName,'WriteMode','append')  
        min_Fu = min(MatrixResults_Fu,[],3);
        writematrix('最小值',saveXlsxName,'WriteMode','append')
        writematrix(min_Fu,saveXlsxName,'WriteMode','append')  
        % % 计算95%置信区间
        % ZX95 = [mean(MatrixResults_Fu,3)-1.96*(std(MatrixResults_Fu,0,3)/sqrt(NumImgs)) mean(MatrixResults_Fu,3)+1.96*(std(MatrixResults_Fu,0,3)/sqrt(NumImgs))];
        % ZX95_Fu = [ZX95(:,1) ZX95(:,6) ZX95(:,2) ZX95(:,7) ZX95(:,3) ZX95(:,8) ZX95(:,4) ZX95(:,9) ZX95(:,5) ZX95(:,10)];
        % writematrix('95%置信区间',saveXlsxName,'WriteMode','append')
        % writematrix(ZX95_Fu,saveXlsxName,'WriteMode','append')  
        

        Matrix_Fu = vertcat(Matrix_Fu,MatrixResults_Fu);
        saveName = fullfile(saveDir,'MatrixAll_Fu.mat'); % saveName = fullfile(saveDir,ErjiDir_list(i_ErjiDir).name,'all.mat');
        save(saveName, 'saveDir','Matrix_Fu');   

    end
    
    fprintf('已保存 %s mat文件！并将该二级目录统计结果打印xlsx \n ', saveName);
end


