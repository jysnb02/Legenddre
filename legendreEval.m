function y = legendreEval(c, x, a, b, normalized)
% y = legendreEval(c, x, a, b, normalized)
% 在点 x 上用系数 c（对应 P_0..P_N）评估多项式
% c: vector length N+1
if nargin < 3 || isempty(a), a = -1; end
if nargin < 4 || isempty(b), b = 1; end
if nargin < 5 || isempty(normalized), normalized = false; end

N = numel(c) - 1;
V = legendreVandermonde(x, N, a, b, normalized);
y = V * c;
end