%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%           LCS_P fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by
%           exploiting local based Pan simulation and golabl injection model
% Interface:
%           I_Fus_LCS_P = LCS_P(I_MS,I_PAN,I_MS_LR,ratio)
%
% Inputs:
%           I_MS:       MS image upsampled at PAN scale;
%           I_PAN:      PAN image;
%           I_MS_LR:    MS image;
%           ratio:      Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           LCS_P:  LCS_P pasharpened image.
%
% References:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_LCS_P = LCS_P(I_MS,I_PAN,ratio)
%% ͨ�������˲����ϳɵͷֱ���ȫɫӰ�� q
r = 4*ratio;%4*ratio;%����Ӱ��ķֱ��ʱ�ֵȷ�����ڿ�ȣ�2r+1
eps = 0;
[q,~,~] = guidedfilter_multibands(I_MS, I_PAN, r, eps);
%% ��ȫɫ�����ȷ�Χ�������ͺϳɵĵͷֱ���Ӱ��һ�£�û��������ΪPan �� q ����ͷǳ��ӽ���
I_PAN = (I_PAN - mean(I_PAN(:)))*std2(q)/std(I_PAN(:)) + mean2(q);
%% ���������˲��ĳɷִ����ں�(CR)��ȫ�ַ���
%choice1��ȫ�ִ���(CR)
I_LCS_P = zeros(size(I_MS));
dtails =(I_PAN-q);
dtails = dtails-mean(dtails(:));
for i = 1:size(I_MS,3)
    %���պϳɵı����ӻ�ȥ
    Gweight = cov(q(:),I_MS(:,:,i))/var(q(:));
    I_LCS_P(:,:,i) = dtails.*Gweight(2) +I_MS(:, :, i);%���Լ�ȥ��ֵ
end
end
