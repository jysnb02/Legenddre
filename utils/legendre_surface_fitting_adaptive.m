function I_correct = legendre_surface_fitting_adaptive(img_gray, param)
    % 勒让德多项式曲面拟合（支持自适应阶数和分块）
    % 
    % 输入参数说明:
    %   img_gray - 输入灰度图像 (double类型)
    %   param - 参数结构体，包含以下字段:
    %       .num_blocks_m - 垂直方向分块数 (推荐: 2-4)
    %       .num_blocks_n - 水平方向分块数 (推荐: 2-4)
    %       .poly_degree - 勒让德多项式阶数 (1-5, 推荐: 3)
    %       .save_blocks - 是否保存中间块图像 (true/false)
    %       .save_path - 块图像保存路径
    %
    % 输出:
    %   I_correct - 拟合后的曲面图像
    
    %% ==================== 参数获取 ====================
    num_blocks_m = param.num_blocks_m;
    num_blocks_n = param.num_blocks_n;
    poly_degree = param.poly_degree;
    
    % 计算每个块的系数数量: (degree+1)*(degree+2)/2
    num_coeffs = (poly_degree + 1) * (poly_degree + 2) / 2;
    total_coeffs = num_coeffs * num_blocks_m * num_blocks_n;
    
    [rows, cols] = size(img_gray);
    
    % 分块边界划分
    split_rows = round(linspace(1, rows, num_blocks_m + 1));
    split_cols = round(linspace(1, cols, num_blocks_n + 1));
    split_rows(end) = rows;
    split_cols(end) = cols;
    
    input_image = double(img_gray);
    [X, Y] = meshgrid(1:cols, 1:rows);
    x_vec = X(:); y_vec = Y(:); z_vec = input_image(:);
    
    %% Step 1: 预计算每个块的基函数和数据
    block_data = cell(num_blocks_m, num_blocks_n);
    G = zeros(total_coeffs, total_coeffs);
    h = zeros(total_coeffs, 1);
    
    for j = 1:num_blocks_m
        for i = 1:num_blocks_n
            x_range = [split_cols(i), split_cols(i+1)];
            y_range = [split_rows(j), split_rows(j+1)];
            
            % 提取块内像素点
            idx = (x_vec >= x_range(1)) & ...
                  ((x_vec < x_range(2)) | ((i == num_blocks_n) & (x_vec == cols))) & ...
                  (y_vec >= y_range(1)) & ...
                  ((y_vec < y_range(2)) | ((j == num_blocks_m) & (y_vec == rows)));
            
            xb = x_vec(idx); yb = y_vec(idx); zb = z_vec(idx);
            
            basis = generate_legendre_basis_adaptive(xb, yb, x_range, y_range, poly_degree);
            B = basis.B;
            
            block_num = (j-1)*num_blocks_n + i;
            coeff_start = (block_num - 1) * num_coeffs;
            
            if ~isempty(B) && ~isempty(zb)
                G_local = B' * B;
                h_local = B' * zb;
                G(coeff_start+1:coeff_start+num_coeffs, coeff_start+1:coeff_start+num_coeffs) = G_local;
                h(coeff_start+1:coeff_start+num_coeffs) = h_local;
            end
            
            block_data{j,i}.x_range = x_range;
            block_data{j,i}.y_range = y_range;
        end
    end
    
    %% Step 2: 构建连续性约束（支持任意分块数）
    constraints = [];
    beq_list = [];
    
    % 垂直边界约束（列间边界）
    if num_blocks_n > 1
        for i = 1:(num_blocks_n - 1)
            x_boundary = split_cols(i+1);
            for j = 1:num_blocks_m
                y_start = split_rows(j);
                y_end = split_rows(j+1);
                if y_end > y_start
                    % 在公共边界上采样
                    yb = (y_start : y_end - 1)';
                    xb = x_boundary * ones(size(yb));
                    
                    left_basis = generate_legendre_basis_adaptive(xb, yb, ...
                        block_data{j,i}.x_range, block_data{j,i}.y_range, poly_degree);
                    right_basis = generate_legendre_basis_adaptive(xb, yb, ...
                        block_data{j,i+1}.x_range, block_data{j,i+1}.y_range, poly_degree);
                    
                    cL_start = ((j-1)*num_blocks_n + (i-1)) * num_coeffs;
                    cR_start = ((j-1)*num_blocks_n + i) * num_coeffs;
                    
                    % 约束数量 = 多项式阶数（G0 到 G_{degree-1}）
                    cons = build_constraints_legendre(left_basis, right_basis, poly_degree, 'x');
                    
                    if ~isempty(cons)
                        A_local = zeros(size(cons,1), total_coeffs);
                        A_local(:, cL_start+1:cL_start+num_coeffs) = cons(:, 1:num_coeffs);
                        A_local(:, cR_start+1:cR_start+num_coeffs) = cons(:, num_coeffs+1:end);
                        constraints = [constraints; A_local];
                        beq_list = [beq_list; zeros(size(cons,1), 1)];
                    end
                end
            end
        end
    end
    
    % 水平边界约束（行间边界）
    if num_blocks_m > 1
        for j = 1:(num_blocks_m - 1)
            y_boundary = split_rows(j+1);
            for i = 1:num_blocks_n
                x_start = split_cols(i);
                x_end = split_cols(i+1);
                if x_end > x_start
                    xb = (x_start : x_end - 1)';
                    yb = y_boundary * ones(size(xb));
                    
                    top_basis = generate_legendre_basis_adaptive(xb, yb, ...
                        block_data{j,i}.x_range, block_data{j,i}.y_range, poly_degree);
                    bottom_basis = generate_legendre_basis_adaptive(xb, yb, ...
                        block_data{j+1,i}.x_range, block_data{j+1,i}.y_range, poly_degree);
                    
                    cT_start = ((j-1)*num_blocks_n + (i-1)) * num_coeffs;
                    cB_start = (j*num_blocks_n + (i-1)) * num_coeffs;
                    
                    cons = build_constraints_legendre(top_basis, bottom_basis, poly_degree, 'y');
                    
                    if ~isempty(cons)
                        A_local = zeros(size(cons,1), total_coeffs);
                        A_local(:, cT_start+1:cT_start+num_coeffs) = cons(:, 1:num_coeffs);
                        A_local(:, cB_start+1:cB_start+num_coeffs) = cons(:, num_coeffs+1:end);
                        constraints = [constraints; A_local];
                        beq_list = [beq_list; zeros(size(cons,1), 1)];
                    end
                end
            end
        end
    end
    
    Aeq = constraints;
    beq = beq_list;
    
    %% Step 3: 求解带约束的最小二乘问题
    if ~isempty(Aeq)
        % 添加小的正则化项提高数值稳定性
        lambda_reg = 1e-6;
        KKT = [G + lambda_reg * speye(size(G)), Aeq'; Aeq, zeros(size(Aeq,1))];
        rhs = [h; beq];
        sol = KKT \ rhs;
        c = sol(1:total_coeffs);
    else
        lambda_reg = 1e-6;
        c = (G + lambda_reg * speye(size(G))) \ h;
    end
    
    %% Step 4: 重建图像
    I_correct = zeros(rows, cols);
    
    do_save = isfield(param, 'save_blocks') && param.save_blocks;
    if do_save
        save_dir = param.save_path;
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
    end
    
    for j = 1:num_blocks_m
        for i = 1:num_blocks_n
            block_num = (j-1)*num_blocks_n + i;
            coeff_range = (block_num - 1) * num_coeffs + (1:num_coeffs);
            
            row_start = split_rows(j);
            if j == num_blocks_m
                row_end = rows;
            else
                row_end = split_rows(j+1) - 1;
            end
            
            col_start = split_cols(i);
            if i == num_blocks_n
                col_end = cols;
            else
                col_end = split_cols(i+1) - 1;
            end
            
            if row_end < row_start || col_end < col_start
                continue;
            end
            
            xb = X(row_start:row_end, col_start:col_end);
            yb = Y(row_start:row_end, col_start:col_end);
            x_range = block_data{j,i}.x_range;
            y_range = block_data{j,i}.y_range;
            
            basis = generate_legendre_basis_adaptive(xb(:), yb(:), x_range, y_range, poly_degree);
            
            if ~isempty(basis.B)
                Z_block = reshape(basis.B * c(coeff_range), size(xb));
                I_correct(row_start:row_end, col_start:col_end) = Z_block;
            else
                I_correct(row_start:row_end, col_start:col_end) = mean(input_image(row_start:row_end, col_start:col_end), 'all');
            end
        end
    end
    
    I_correct(isnan(I_correct)) = mean(input_image(:));
end