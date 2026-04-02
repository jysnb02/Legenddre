function S = L0Smoothing(I, lambda, beta_max, kappa, iter_max)
% 参数设置
if nargin < 2, lambda = 0.02; end
if nargin < 3, beta_max = 1e5; end
if nargin < 4, kappa = 2.0; end
if nargin < 5, iter_max = 20; end

% 转换灰度图并归一化
if size(I, 3) == 3, I = rgb2gray(I); end
I = im2double(I);
[m, n] = size(I);

% 初始化变量
S = I;
beta = 2*lambda;
h = zeros(m, n);
v = zeros(m, n);

% 迭代优化
for iter = 1:iter_max
    % 计算梯度项的转置卷积
    DxTh = imfilter(h, [-1, 1], 'symmetric'); % 水平方向后向差分
    DyTv = imfilter(v, [-1; 1], 'symmetric'); % 垂直方向后向差分
    
    % 构造线性方程右边项
    rhs = I + beta*(DxTh + DyTv);
    
    % DCT变换求解
    dct_rhs = dct2(rhs);
    [x, y] = meshgrid(0:n-1, 0:m-1);
    den = 1 + beta*(2 - 2*cos(pi*y/m) + 2 - 2*cos(pi*x/n));
    den(1,1) = 1; % 避免除以0
    
    % 更新平滑图像
    S = idct2(dct_rhs ./ den);
    
    % 计算当前梯度
    dx = imfilter(S, [1, -1], 'symmetric'); % 水平前向差分
    dy = imfilter(S, [1; -1], 'symmetric'); % 垂直前向差分
    
    % 更新辅助变量
    threshold = lambda / beta;
    mask = (dx.^2 + dy.^2) >= threshold;
    h = mask .* dx;
    v = mask .* dy;
    
    % 更新beta参数
    beta = min(beta*kappa, beta_max);
end

% 使用示例
% img = imread('input.jpg');
% S = L0Smoothing(img);
% imshow([img, S]);