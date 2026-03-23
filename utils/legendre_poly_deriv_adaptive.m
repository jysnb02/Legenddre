function L_deriv = legendre_poly_deriv_adaptive(degree, x, order)
    % 计算勒让德多项式的 order 阶导数
    % 输入:
    %   degree - 最高阶数
    %   x - 输入向量
    %   order - 导数阶数 (0,1,2,3,4,5,...)
    % 输出:
    %   L_deriv - (n x (degree+1)) 矩阵
    
    if order == 0
        L_deriv = legendre_poly_adaptive(degree, x);
        return;
    end
    
    n = length(x);
    L_deriv = zeros(n, degree + 1);
    
    if order == 1
        % 一阶导数递推公式: L'_0=0, L'_1=1, L'_k = (2k-1)L_{k-1} + L'_{k-2}
        if degree >= 1
            L_deriv(:, 2) = 1;
            L = legendre_poly_adaptive(degree, x);
            for k = 2:degree
                L_deriv(:, k+1) = (2*k - 1) * L(:, k) + L_deriv(:, k);
            end
        end
        
    elseif order == 2
        % 二阶导数: 使用差分法计算
        if degree >= 2
            eps = 1e-8;
            for k = 2:degree
                x_plus = x + eps;
                x_minus = x - eps;
                L1_plus = legendre_poly_deriv_adaptive(k, x_plus, 1);
                L1_minus = legendre_poly_deriv_adaptive(k, x_minus, 1);
                L_deriv(:, k+1) = (L1_plus(:, k+1) - L1_minus(:, k+1)) / (2 * eps);
            end
        end
        
    elseif order == 3
        % 三阶导数
        if degree >= 3
            eps = 1e-8;
            for k = 3:degree
                x_plus = x + eps;
                x_minus = x - eps;
                L2_plus = legendre_poly_deriv_adaptive(k, x_plus, 2);
                L2_minus = legendre_poly_deriv_adaptive(k, x_minus, 2);
                L_deriv(:, k+1) = (L2_plus(:, k+1) - L2_minus(:, k+1)) / (2 * eps);
            end
        end
        
    elseif order == 4
        % 四阶导数
        if degree >= 4
            eps = 1e-8;
            for k = 4:degree
                x_plus = x + eps;
                x_minus = x - eps;
                L3_plus = legendre_poly_deriv_adaptive(k, x_plus, 3);
                L3_minus = legendre_poly_deriv_adaptive(k, x_minus, 3);
                L_deriv(:, k+1) = (L3_plus(:, k+1) - L3_minus(:, k+1)) / (2 * eps);
            end
        end
        
    elseif order == 5
        % 五阶导数
        if degree >= 5
            eps = 1e-8;
            for k = 5:degree
                x_plus = x + eps;
                x_minus = x - eps;
                L4_plus = legendre_poly_deriv_adaptive(k, x_plus, 4);
                L4_minus = legendre_poly_deriv_adaptive(k, x_minus, 4);
                L_deriv(:, k+1) = (L4_plus(:, k+1) - L4_minus(:, k+1)) / (2 * eps);
            end
        end
    end
end