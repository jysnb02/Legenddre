function c = legendreFit(x, y, N, a, b, normalized, weights)
% c = legendreFit(x, y, N, a, b, normalized, weights)
% 用勒让德多项式拟合数据 (x,y)，返回系数 c (length N+1)，使得
% y ≈ V * c，其中 V 列为 P_0..P_N
% weights: 可选，向量或对角权重（加权最小二乘）
% 
if nargin < 4, a = -1; end
if nargin < 5, b = 1; end
if nargin < 6 || isempty(normalized), normalized = false; end
if nargin < 7, weights = []; end

V = legendreVandermonde(x, N, a, b, normalized);

if isempty(weights)
    % 普通最小二乘
    c = V \ y;
else
    w = weights(:);
    if numel(w) ~= size(V,1)
        error('weights must have same length as x');
    end
    % solve (W*V) c = W*y  where W = diag(sqrt(w))
    W = diag(sqrt(w));
    c = (W*V) \ (W*y);
end
end