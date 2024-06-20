function varargout = guidedfilter_multibands(I, p, r, eps)
%     guidedfilter_multibands O(1) time implementation of guided filter using a color image as the guidance.
%   - guidance image: I (should be a image with multi-bands (gt 1))
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps
% Create by Leiguang wang at wlgbain@126.com in 2019/03/10 based on GUIDEDFILTER_COLOR
%

[hei, wid,dim] = size(I);
% [hei, wid] = size(p);

N = boxfilter(ones(hei, wid), r); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.

mean_I = zeros([hei, wid,dim]);
for i = 1:dim
    mean_I(:,:,i) = boxfilter(I(:, :, i), r) ./ N;
end

mean_p = boxfilter(p, r) ./ N;

mean_Ip = zeros([hei, wid,dim]);
for i = 1:dim
    mean_Ip(:,:,i) = boxfilter(I(:, :, i).*p, r) ./ N;
end

% covariance of (I, p) in each local patch.
cov_Ip = zeros([hei, wid,dim]);
for i = 1:dim
    cov_Ip(:,:,i) =  mean_Ip(:,:,i)-mean_I(:,:,i).*mean_p;
end

% variance of I in each local patch: the matrix Sigma in Eqn (14).
% Note the variance in each local patch is a bandNum x bandNum symmetric matrix:
%           rr, rg, rb
%   Sigma = rg, gg, gb
%           rb, gb, bb
var_I= zeros(hei, wid,dim*dim);
Idx = 1;
for i = 1:dim
    for j = 1:dim
        var_I(:,:,Idx) =boxfilter(I(:, :, i).*I(:, :, j), r) ./ N - mean_I(:,:,i) .*  mean_I(:,:,j);
        Idx = Idx+1;
    end
end

%{
%M2 ���м��㣬��ʱ�Գ�,ռ���ڴ��
a = zeros(hei*wid, dim);
tmp1 = reshape(permute(var_I,[3,1,2]),[dim,dim,hei*wid])+repmat(eps * eye(dim),[1,1,hei*wid]);
tmp2 = reshape(permute(cov_Ip,[3,1,2]),[1,dim,hei*wid]);
parfor i = 1:hei*wid
%  for i = 1:hei*wid
    a(i,:) = tmp2(:,:,i)/tmp1(:,:,i);% Eqn. (14) in the paper;cov_Ip_Local * inv(Sigma + eps * eye(dim));
end
a = reshape(a,[hei,wid, dim]);

%M3 ����mex,��ʱ����,ռ���ڴ��
var_I = reshape(var_I,[hei*wid,dim*dim]);
cov_Ip = reshape(cov_Ip,[hei*wid,dim]);
a = CalcuateA( var_I,cov_Ip, hei,wid, dim,eps);
a = permute(reshape(a,[hei,wid, dim]),[2,1,3]);

%M1��ʱ��,ռ���ڴ�С
a = zeros(hei, wid, dim);
for y=1:hei
    for x=1:wid
        Sigma = reshape(var_I(y,x,:), [dim,dim]);
       cov_Ip_Local = reshape(cov_Ip(y, x,:), [1,dim]);
        a(y, x, :) = cov_Ip_Local/(Sigma + eps * eye(dim)); % Eqn. (14) in the paper;cov_Ip_Local * inv(Sigma + eps * eye(dim));
    end
end

%M4 ����mex,��ʱ��Ȼ�ܳ�,ռ���ڴ�Ҳ�������⣡��
var_I = reshape(var_I,[hei*wid,dim*dim]);
cov_Ip = reshape(cov_Ip,[hei*wid,dim]);
a = CalcuateB(var_I,cov_Ip, hei,wid, dim,eps);
a = reshape(a,[hei,wid, dim]);

%M5 :�ֿ鴦����ʱ����ƫ��
tmp1 = reshape(permute(var_I,[3,1,2]),[dim,dim,hei*wid])+repmat(eps * eye(dim),[1,1,hei*wid]);
tmp2 = reshape(permute(cov_Ip,[3,1,2]),[1,dim,hei*wid]);
bs = cat(1,tmp1,tmp2);
bs = reshape(bs, [dim+1,dim*hei*wid]);
fun_eg = @(bs) calcuLocA(bs.data, dim);
a = blockproc(bs,[dim+1 dim],fun_eg);
a = reshape(a, [dim,wid,hei]);
a = permute(a,[2,3,1]);

%M6:��Ӱ�����������Э�������ֲ�Э����, ��΢���;��ȵķ�������߼���Ч�ʣ�Լ5���������Ǿ��Ƚ���
cov_Ip = reshape(cov_Ip,[hei*wid,dim]);
a = cov_Ip/(cov(reshape(I,[hei*wid, dim])));% Eqn. (14) in the paper;cov_Ip_Local * inv(Sigma + eps * eye(dim));% ����ֲ�Э�������MS����Э����,��������
a = permute(reshape(a,[hei,wid, dim]),[2,1,3]);

%}

%����СӰ�񣬲���parfor Ч����ߣ�for ��֮����Ӱ��forЧ�����
a = zeros(hei*wid, dim);
tmp1 = reshape(permute(var_I,[3,1,2]),[dim,dim,hei*wid])+repmat(eps * eye(dim),[1,1,hei*wid]);%��������
tmp2 = reshape(permute(cov_Ip,[3,1,2]),[1,dim,hei*wid]);
% parfor i = 1:hei*wid
% index =[];
for i = 1:hei*wid
%     if cond(tmp1(:,:,i)) <100000
        a(i,:) = tmp2(:,:,i)/tmp1(:,:,i);% Eqn. (14) in the paper;cov_Ip_Local * inv(Sigma + eps * eye(dim));
%     else
%         index = cat(1,index,i);
%     end
end
a = reshape(a,[hei,wid, dim]);

q = 0;
for i = 1:dim
    %������ܳ��ֵ�tmp2����������������������Ȩ��(��Ҫ���QB_SWFU)
    Flag = find(isnan(a(:,:,i)));
    if ~isempty(Flag)
        tmp_a = a(:,:,i);
        tmp_I = I(:,:,i)./sum(I,3);
        tmp_a(Flag) = tmp_I(Flag);
        a(:,:,i) = tmp_a;
    end
        
    a(:, :, i) = boxfilter(a(:, :, i), r)./N;%
    q = q+a(:, :, i).* I(:, :, i);
end

%
b = mean_p - sum(a.*mean_I,3);
b = boxfilter(b, r)./N; %ȡƽ��

%�˲����
q = q+b;

if nargout == 1
    varargout{1} = q;
elseif nargout == 2
    varargout{1} = q;
    varargout{2} = a;
elseif nargout == 3
    varargout{1} = q;
    varargout{2} = a;
    varargout{3} = b;
end

function a = calcuLocA(bs, dims)
    a = (bs(end,:)/bs(1:dims,:));
return