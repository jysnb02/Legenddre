function basis = generate_legendre_basis_adaptive(x, y, x_range, y_range, degree)
    % 生成勒让德基函数及其导数（支持任意阶数）
    % 
    % 输入:
    %   x, y - 像素坐标向量
    %   x_range, y_range - 块的坐标范围
    %   degree - 多项式阶数
    % 输出:
    %   basis - 包含基函数和导数的结构体
    
    num_pts = length(x);
    if num_pts == 0
        basis = struct('B', [], 'Bx', [], 'By', [], 'Bxx', [], 'Byy', [], 'Bxy', [], ...
                      'Bxxx', [], 'Byyy', [], 'Bxxy', [], 'Bxyy', []);
        return;
    end
    
    % 坐标归一化到 [-1, 1]
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
    
    % 计算勒让德多项式及其导数
    Lx = legendre_poly_adaptive(degree, x_scaled);
    Ly = legendre_poly_adaptive(degree, y_scaled);
    Lx_d1 = legendre_poly_deriv_adaptive(degree, x_scaled, 1);
    Ly_d1 = legendre_poly_deriv_adaptive(degree, y_scaled, 1);
    Lx_d2 = legendre_poly_deriv_adaptive(degree, x_scaled, 2);
    Ly_d2 = legendre_poly_deriv_adaptive(degree, y_scaled, 2);
    Lx_d3 = legendre_poly_deriv_adaptive(degree, x_scaled, 3);
    Ly_d3 = legendre_poly_deriv_adaptive(degree, y_scaled, 3);
    
    % 系数数量: 对于阶数 degree，系数个数 = (degree+1)*(degree+2)/2
    num_coeffs = (degree + 1) * (degree + 2) / 2;
    
    % 初始化基函数矩阵
    B = zeros(num_pts, num_coeffs);
    Bx = zeros(num_pts, num_coeffs);
    By = zeros(num_pts, num_coeffs);
    Bxx = zeros(num_pts, num_coeffs);
    Byy = zeros(num_pts, num_coeffs);
    Bxy = zeros(num_pts, num_coeffs);
    Bxxx = zeros(num_pts, num_coeffs);
    Byyy = zeros(num_pts, num_coeffs);
    Bxxy = zeros(num_pts, num_coeffs);
    Bxyy = zeros(num_pts, num_coeffs);
    
    % 构造基函数 (按顺序: i=0..degree, j=0..degree-i)
    k = 0;
    for i = 0:degree
        for j = 0:(degree - i)
            k = k + 1;
            
            % 函数值基
            B(:, k) = Lx(:, i+1) .* Ly(:, j+1);
            
            % 一阶导数基
            Bx(:, k) = (Lx_d1(:, i+1) .* Ly(:, j+1)) * dx_scaling;
            By(:, k) = (Lx(:, i+1) .* Ly_d1(:, j+1)) * dy_scaling;
            
            % 二阶导数基
            Bxx(:, k) = (Lx_d2(:, i+1) .* Ly(:, j+1)) * dx_scaling^2;
            Byy(:, k) = (Lx(:, i+1) .* Ly_d2(:, j+1)) * dy_scaling^2;
            Bxy(:, k) = (Lx_d1(:, i+1) .* Ly_d1(:, j+1)) * dx_scaling * dy_scaling;
            
            % 三阶导数基（用于高阶约束）
            Bxxx(:, k) = (Lx_d3(:, i+1) .* Ly(:, j+1)) * dx_scaling^3;
            Byyy(:, k) = (Lx(:, i+1) .* Ly_d3(:, j+1)) * dy_scaling^3;
            Bxxy(:, k) = (Lx_d2(:, i+1) .* Ly_d1(:, j+1)) * dx_scaling^2 * dy_scaling;
            Bxyy(:, k) = (Lx_d1(:, i+1) .* Ly_d2(:, j+1)) * dx_scaling * dy_scaling^2;
        end
    end
    
    basis = struct('B', B, 'Bx', Bx, 'By', By, ...
                   'Bxx', Bxx, 'Byy', Byy, 'Bxy', Bxy, ...
                   'Bxxx', Bxxx, 'Byyy', Byyy, 'Bxxy', Bxxy, 'Bxyy', Bxyy);
end