function V = legendreVandermonde(x, N, a, b, normalized)
% V = legendreVandermonde(x, N, a, b, normalized)
% 构造勒让德多项式的 Vandermonde 矩阵，列为 P_0 ... P_N（degree 0..N）
% x: 向量
% N: 最大多项式次数（非项数）
% a,b: 可选，数据原始区间，默认 [-1,1]。若给出，���先把 x 映射到 [-1,1]。
% normalized: 可选布尔，若 true 则对每列做正交归一化（orthonormal）。
%
% 使用递推：(n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}

if nargin < 3 || isempty(a), a = -1; end
if nargin < 4 || isempty(b), b = 1; end
if nargin < 5 || isempty(normalized), normalized = false; end

x = x(:);
m = numel(x);

% 映射到 [-1,1]
if ~(a == -1 && b == 1)
    xs = (2*(x - a)/(b - a) - 1);
else
    xs = x;
end

V = zeros(m, N+1);
P0 = ones(m,1);
V(:,1) = P0;
if N >= 1
    P1 = xs;
    V(:,2) = P1;
end

Pnm1 = P0; % P_{n-1}
Pn = [];
if N >= 1
    Pn = P1; % P_n
end

for n = 1:(N-1)
    % compute P_{n+1}
    % (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
    Pnp1 = ((2*n+1) .* xs .* Pn - n .* Pnm1) ./ (n+1);
    V(:, n+2) = Pnp1;
    Pnm1 = Pn;
    Pn = Pnp1;
end

if normalized
    % 标准勒让德正交关系： <P_n, P_m> = 2/(2n+1) delta_{nm}
    % 要变成正交归一化 multiply each P_n by sqrt((2n+1)/2)
    for n = 0:N
        scale = sqrt((2*n+1)/2);
        V(:, n+1) = V(:, n+1) * scale;
    end
end
end
