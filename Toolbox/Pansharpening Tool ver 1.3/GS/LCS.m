%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           LCS fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting local based Pan simulation and lcoal injection
% 
% Interface:
%           I_Fus_LCS_P = LCS(I_MS,I_PAN,I_MS_LR,ratio)
%
% Inputs:
%           I_MS:       MS image upsampled at PAN scale;
%           I_PAN:      PAN image;
%           ratio:      Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           LCS:  LCS pasharpened image.
% 
% References:
%           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_LCS = LCS(I_MS,I_PAN,ratio)
%% ͨ�������˲����ϳɵͷֱ���ȫɫӰ�� q
r = 4*ratio;%����Ӱ��ķֱ��ʱ�ֵȷ�����ڿ�ȣ�2r+1
eps = 0;
[q,~,~] = guidedfilter_multibands(I_MS, I_PAN, r, eps);
%% ��ȫɫ�����ȷ�Χ�������ͺϳɵĵͷֱ���Ӱ��һ�£�û��������ΪPan �� q ����ͷǳ��ӽ���
I_PAN = (I_PAN - mean(I_PAN(:)))*std2(q)/std(I_PAN(:)) + mean2(q);
%% ���������˲��ĳɷִ����ں�(CR)���ֲ�����
dtails =(I_PAN-q);
% dtails = dtails-mean(dtails(:));%�Ƿ��ȥ��ֵ�������ս��Ӱ�첻��

I_LCS = zeros(size(I_MS));
[hei, wid, ~]= size(I_MS);
N = boxfilter(ones(hei, wid), r); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.
for i = 1:size(I_MS,3)    
        %     %M5:�ֲ�GS����:�������һ������ȥ��ֵ���ټ��뵽ԭʼӰ����
    weight = calLcoalCova(q,I_MS(:,:,i),r)./(calLcoalCova(q,q,r));
    weight = boxfilter(weight, r)./N;
    I_LCS(:,:,i) = dtails.*weight +I_MS(:, :, i);  
    %�����ò����������ںϾ��ȵĸ���
%      I_LCS(:,:,i) = I_LCS(:,:,i)-mean2(I_LCS(:,:,i))+mean2(I_MS(:, :,i));
%      
end
end
