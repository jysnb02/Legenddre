function L = legendre_poly_adaptive(degree, x)
    % 计算 0 到 degree 阶勒让德多项式在 x 处的值
    % 递推公式: L0=1, L1=x, (k+1)L_{k+1} = (2k+1)x L_k - k L_{k-1}
    %
    % 输入:
    %   degree - 最高阶数
    %   x - 输入向量 (已归一化到 [-1,1])
    % 输出:
    %   L - (n x (degree+1)) 矩阵，第 k 列为 k-1 阶多项式值
    
    n = length(x);
    L = zeros(n, degree + 1);
    
    L(:, 1) = 1;  % L0(x) = 1
    
    if degree >= 1
        L(:, 2) = x;  % L1(x) = x
        for k = 2:degree
            % 递推公式: (k+1)L_{k+1} = (2k+1)x L_k - k L_{k-1}
            L(:, k+1) = ((2*k - 1) * x .* L(:, k) - (k - 1) * L(:, k-1)) / k;
        end
    end
end