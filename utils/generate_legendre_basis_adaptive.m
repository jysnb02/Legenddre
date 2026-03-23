function basis = generate_legendre_basis_adaptive(x, y, x_range, y_range, degree)
    % 生成自适应阶数的勒让德基函数
    % 输入:
    %   x, y - 像素坐标
    %   x_range, y_range - 块坐标范围
    %   degree - 多项式阶数
    % 输出:
    %   basis - 包含基函数和导数的结构体
    
    num_pts = length(x);
    if num_pts == 0
        basis = struct('B', []);
        return;
    end
    
    % 坐标归一化
    if x_range(2) > x_range(1)
        x_scaled = 2 * (x - x_range(1)) / (x_range(2) - x_range(1)) - 1;
        dx_scaling = 2 / (x_range(2) - x_range(1));
    else
        x_scaled = zeros(num_pts, 1);
        dx_scaling = 0;
    end
    
    if y_range(2) > y_range(1)
        y_scaled = 2 * (y - y_range(1)) / (y_range(2) - y_range(1)) - 1;
        dy_scaling = 2 / (y_range(2) - y_range(1));
    else
        y_scaled = zeros(num_pts, 1);
        dy_scaling = 0;
    end
    
    % 预计算所有导数 (0 到 degree-1 阶)
    num_derivs = degree;
    Lx_cell = cell(1, num_derivs + 1);
    Ly_cell = cell(1, num_derivs + 1);
    
    for d = 0:num_derivs
        Lx_cell{d+1} = legendre_poly_deriv_adaptive(degree, x_scaled, d);
        Ly_cell{d+1} = legendre_poly_deriv_adaptive(degree, y_scaled, d);
    end
    
    % 系数数量
    num_coeffs = (degree + 1) * (degree + 2) / 2;
    
    % 初始化结构体
    basis = struct();
    basis.B = zeros(num_pts, num_coeffs);
    
    % 创建导数场
    deriv_names_x = {'Bx', 'Bxx', 'Bxxx', 'Bxxxx', 'Bxxxxx'};
    deriv_names_y = {'By', 'Byy', 'Byyy', 'Byyyy', 'Byyyyy'};
    
    for d = 1:degree
        if d <= length(deriv_names_x)
            basis.(deriv_names_x{d}) = zeros(num_pts, num_coeffs);
            basis.(deriv_names_y{d}) = zeros(num_pts, num_coeffs);
        end
    end
    
    % 构造基函数
    k = 0;
    for i = 0:degree
        for j = 0:(degree - i)
            k = k + 1;
            
            % 基础基函数
            basis.B(:, k) = Lx_cell{1}(:, i+1) .* Ly_cell{1}(:, j+1);
            
            % x方向导数
            for d = 1:degree
                if d <= i && d <= length(deriv_names_x)
                    basis.(deriv_names_x{d})(:, k) = (Lx_cell{d+1}(:, i+1) .* Ly_cell{1}(:, j+1)) * (dx_scaling^d);
                end
            end
            
            % y方向导数
            for d = 1:degree
                if d <= j && d <= length(deriv_names_y)
                    basis.(deriv_names_y{d})(:, k) = (Lx_cell{1}(:, i+1) .* Ly_cell{d+1}(:, j+1)) * (dy_scaling^d);
                end
            end
        end
    end
end