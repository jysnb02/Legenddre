function L = legendre_poly_adaptive(degree, x)
    % 计算 0 到 degree 阶勒让德多项式
    % 输入:
    %   degree - 最高阶数
    %   x - 输入向量 (已归一化到 [-1,1])
    % 输出:
    %   L - (n x (degree+1)) 矩阵
    
    n = length(x);
    L = zeros(n, degree + 1);
    
    L(:, 1) = 1;  % L0(x)=1
    
    if degree >= 1
        L(:, 2) = x;  % L1(x)=x
        for k = 2:degree
            % 递推公式: (k+1)L_{k+1}(x) = (2k+1)xL_k(x) - kL_{k-1}(x)
            L(:, k+1) = ((2*k - 1) * x .* L(:, k) - (k - 1) * L(:, k-1)) / k;
        end
    end
end