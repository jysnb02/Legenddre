function cons = build_constraints_adaptive(basis1, basis2, degree, direction)
    % 根据阶数动态构造连续性约束
    % 输入:
    %   basis1, basis2 - 相邻块的基函数结构体
    %   degree - 多项式阶数
    %   direction - 方向: 'x' 或 'y'
    % 输出:
    %   cons - 约束矩阵
    
    num_coeffs = size(basis1.B, 2);
    num_constraints = degree;  % G0 到 G_{degree-1}
    num_pts = size(basis1.B, 1);
    
    cons = zeros(num_constraints * num_pts, 2 * num_coeffs);
    
    % 根据方向选择导数字段
    if strcmp(direction, 'x')
        deriv_fields = {'B', 'Bx', 'Bxx', 'Bxxx', 'Bxxxx', 'Bxxxxx'};
    else
        deriv_fields = {'B', 'By', 'Byy', 'Byyy', 'Byyyy', 'Byyyyy'};
    end
    
    row_offset = 0;
    for d = 1:degree
        field_name = deriv_fields{d};
        
        % 获取对应导数的基函数
        if isfield(basis1, field_name) && isfield(basis2, field_name)
            B1 = basis1.(field_name);
            B2 = basis2.(field_name);
        else
            cons = [];
            return;
        end
        
        % 填充约束矩阵
        cons(row_offset+1:row_offset+num_pts, 1:num_coeffs) = B1;
        cons(row_offset+1:row_offset+num_pts, num_coeffs+1:2*num_coeffs) = -B2;
        
        row_offset = row_offset + num_pts;
    end
end