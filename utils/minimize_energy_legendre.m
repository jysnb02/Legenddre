function [fb, b, cost_history] = minimize_energy_legendre(g, varargin)
    % 能量最小化主函数（勒让德版本）
    % min_{b,fb} ||g - fb - b||_2^2 + η||Dg - Dfb - Db||_2^2
    
    %% 参数设置
    params = parse_inputs_simple(varargin{:});
    
    defaults.eta = 0.1;
    defaults.max_iter = 50;
    defaults.tol = 1e-6;
    defaults.display = true;
    defaults.poly_degree = 3;
    defaults.num_blocks_m = 2;
    defaults.num_blocks_n = 2;
    
    fields = fieldnames(defaults);
    for i = 1:length(fields)
        if ~isfield(params, fields{i})
            params.(fields{i}) = defaults.(fields{i});
        end
    end
    
    %% 初始化
    [M, N] = size(g);
    
    param_fit.num_blocks_m = params.num_blocks_m;
    param_fit.num_blocks_n = params.num_blocks_n;
    param_fit.poly_degree = params.poly_degree;
    param_fit.save_blocks = false;
    
    % 使用勒让德曲面拟合初始化
    b_initial = legendre_surface_fitting_adaptive(g, param_fit);
    
    b = b_initial;
    fb = g - b;
    
    % 梯度算子
    [Dx, Dy] = create_gradient_operators_centered(M, N);
    [Dgx, Dgy] = gradient(g);
    Dg = Dgx + 1i * Dgy;
    
    cost_history = zeros(params.max_iter, 1);
    
    %% 交替优化
    for iter = 1:params.max_iter
        fb_prev = fb;
        b_prev = b;
        
        % 固定b，优化fb
        [A_fb, rhs_fb] = build_fb_system_simple(g, b, Dg, Dx, Dy, params.eta, M, N);
        A_fb_reg = A_fb + 1e-8 * speye(size(A_fb));
        fb_vec = A_fb_reg \ rhs_fb(:);
        fb = reshape(fb_vec, M, N);
        
        % 固定fb，优化b（使用勒让德拟合）
        [A_b, rhs_b] = build_b_system_simple(g, fb, Dg, Dx, Dy, params.eta, M, N);
        A_b_reg = A_b + 1e-8 * speye(size(A_b));
        b_vec = A_b_reg \ rhs_b(:);
        b = reshape(b_vec, M, N);
        
        % 勒让德平滑
        b = legendre_surface_fitting_adaptive(b, param_fit);
        
        % 成本函数
        cost = compute_cost_simple(g, fb, b, Dg, Dx, Dy, params.eta);
        cost_history(iter) = cost;
        
        % 收敛检测
        fb_change = norm(fb(:) - fb_prev(:)) / (norm(fb_prev(:)) + eps);
        b_change = norm(b(:) - b_prev(:)) / (norm(b_prev(:)) + eps);
        
        if params.display
            fprintf('迭代 %3d: 成本 = %.6e, fb变化 = %.6e, b变化 = %.6e\n', ...
                    iter, cost, fb_change, b_change);
        end
        
        if max(fb_change, b_change) < params.tol
            if params.display
                fprintf('在迭代 %d 收敛\n', iter);
            end
            cost_history = cost_history(1:iter);
            break;
        end
    end
    
    if params.display
        visualize_results_simple(g, fb, b, cost_history);
    end
end

function [Dx, Dy] = create_gradient_operators_centered(M, N)
    total_pixels = M * N;
    
    % 水平梯度
    i_x = []; j_x = []; s_x = [];
    for row = 1:M
        for col = 2:(N-1)
            idx = (row-1)*N + col;
            i_x = [i_x; idx; idx];
            j_x = [j_x; idx-1; idx+1];
            s_x = [s_x; -0.5; 0.5];
        end
    end
    Dx = sparse(i_x, j_x, s_x, total_pixels, total_pixels);
    
    % 垂直梯度
    i_y = []; j_y = []; s_y = [];
    for row = 2:(M-1)
        for col = 1:N
            idx = (row-1)*N + col;
            i_y = [i_y; idx; idx];
            j_y = [j_y; idx-N; idx+N];
            s_y = [s_y; -0.5; 0.5];
        end
    end
    Dy = sparse(i_y, j_y, s_y, total_pixels, total_pixels);
end

function [A, rhs] = build_fb_system_simple(g, b, Dg, Dx, Dy, eta, M, N)
    total_pixels = M * N;
    I = speye(total_pixels);
    DTD = Dx' * Dx + Dy' * Dy;
    A = I + eta * DTD;
    rhs = (g(:) - b(:)) + eta * (Dx' * real(Dg(:) - (Dx*b(:) + 1i*Dy*b(:))) + ...
                                   Dy' * imag(Dg(:) - (Dx*b(:) + 1i*Dy*b(:))));
end

function [A, rhs] = build_b_system_simple(g, fb, Dg, Dx, Dy, eta, M, N)
    total_pixels = M * N;
    I = speye(total_pixels);
    DTD = Dx' * Dx + Dy' * Dy;
    A = I + eta * DTD;
    rhs = (g(:) - fb(:)) + eta * (Dx' * real(Dg(:) - (Dx*fb(:) + 1i*Dy*fb(:))) + ...
                                   Dy' * imag(Dg(:) - (Dx*fb(:) + 1i*Dy*fb(:))));
end

function params = parse_inputs_simple(varargin)
    params = struct;
    for i = 1:2:length(varargin)
        if i+1 <= length(varargin)
            params.(varargin{i}) = varargin{i+1};
        end
    end
end

function cost = compute_cost_simple(g, fb, b, Dg, Dx, Dy, eta)
    [M, N] = size(g);
    data_fidelity = sum((g(:) - fb(:) - b(:)).^2);
    Dfb = Dx*fb(:) + 1i*Dy*fb(:);
    Db = Dx*b(:) + 1i*Dy*b(:);
    grad_fidelity = eta * sum(abs(Dg(:) - Dfb - Db).^2);
    cost = data_fidelity + grad_fidelity;
end

function visualize_results_simple(g, fb, b, cost_history)
    figure('Position', [150, 150, 1200, 600]);
    
    subplot(2,3,1); imshow(g, []); title('原始图像'); colorbar;
    subplot(2,3,2); imshow(fb, []); title('目标场景'); colorbar;
    subplot(2,3,3); imshow(b, []); title('热辐射效应'); colorbar;
    subplot(2,3,6); surf(fliplr(b)); title('热辐射效应 3D'); shading interp; axis tight; axis off;
    subplot(2,3,4); imshow(fb+b, []); title('重建结果'); colorbar;
    subplot(2,3,5); semilogy(1:length(cost_history), cost_history, 'LineWidth', 2);
    xlabel('迭代次数'); ylabel('成本函数值'); title('收敛曲线'); grid on;
end