%% ============================================================
% 勒让德多项式热辐射偏置场校正主程序
% 功能：批量处理图像，进行勒让德曲面拟合，去除热辐射偏置场
% 支持分块数量和多项式阶数可调
% ============================================================

clc;
clear;
close all;
addpath('./utils/');

%% ==================== 参数配置区域（可调参数） ====================
% 用户可根据需要修改以下参数

% ----- 1. 文件路径配置 -----
config.input_folder = 'E:\test_legenddre\IRBFD\退化图';      % 输入图像文件夹
config.output_folder = 'E:\test_legenddre\res\IRBFD';  % 输出文件夹
config.clearer_folder = 'E:\test_legenddre\IRBFD\清晰图';    % 参考图像文件夹（可选）

% ----- 2. 勒让德拟合参数 -----
legendre_param.num_blocks_m = 3;      % 垂直方向分块数（推荐：2-4）
legendre_param.num_blocks_n = 3;      % 水平方向分块数（推荐：2-4）
legendre_param.poly_degree = 4;       % 勒让德多项式阶数（1-5，推荐：3）
legendre_param.save_blocks = true;    % 是否保存中间块图像
legendre_param.save_path = config.output_folder;  % 块图像保存路径

% ----- 3. 能量最小化参数 -----
energy_param.eta = 0.001;              % 梯度保真项权重（推荐：0.001-0.1）
energy_param.max_iter = 5;            % 最大迭代次数（推荐：3-10）
energy_param.tol = 1e-6;              % 收敛容差
energy_param.display = true;          % 是否显示迭代信息

% ----- 4. 处理选项 -----
use_energy_minimization = false;       % true: 交替优化, false: 直接拟合
save_results = true;                  % 是否保存结果
show_results = true;                  % 是否显示结果图

%% ==================== 参数说明（供用户参考） ====================
fprintf('\n');
fprintf('==================== 参数说明 ====================\n');
fprintf('1. 分块数 (num_blocks_m, num_blocks_n):\n');
fprintf('   - 越大越能捕捉局部细节，但计算量增加\n');
fprintf('   - 推荐: 2×2 或 3×3\n\n');
fprintf('2. 多项式阶数 (poly_degree):\n');
fprintf('   - 1阶: 平面拟合（最简单）\n');
fprintf('   - 2阶: 二次曲面（适合平滑变化）\n');
fprintf('   - 3阶: 三次曲面（推荐，平衡精度和稳定性）\n');
fprintf('   - 4阶: 四次曲面（更精细，可能过拟合）\n');
fprintf('   - 5阶: 五次曲面（需要更多数据点）\n\n');
fprintf('3. 梯度权重 (eta):\n');
fprintf('   - 值越大，梯度约束越强，图像越平滑\n');
fprintf('   - 推荐: 0.001-0.1\n\n');
fprintf('4. 迭代次数 (max_iter):\n');
fprintf('   - 越多效果越好，但耗时增加\n');
fprintf('   - 推荐: 3-10\n');
fprintf('====================================================\n\n');

%% ==================== 创建输出文件夹 ====================
if ~exist(config.output_folder, 'dir')
    mkdir(config.output_folder);
end

%% ==================== 获取输入图像列表 ====================
file_list = dir(fullfile(config.input_folder, '*.png'));
if isempty(file_list)
    file_list = dir(fullfile(config.input_folder, '*.jpg'));
end
if isempty(file_list)
    file_list = dir(fullfile(config.input_folder, '*.bmp'));
end
if isempty(file_list)
    file_list = dir(fullfile(config.input_folder, '*.tif'));
end

if isempty(file_list)
    error('未找到图像文件，请检查路径: %s', config.input_folder);
end

fprintf('找到 %d 张图像\n\n', length(file_list));

%% ==================== 初始化结果存储 ====================
results = cell(length(file_list), 12);
result_names = {'图像名', 'PSNR', 'SSIM', 'EOG', 'NIQE', 'NISA', ...
                'NIQE_orig', 'NISA_orig', 'PSNR_ref', 'SSIM_ref', ...
                '耗时(s)', '阶数'};

%% ==================== 主循环：处理每一张图像 ====================
for idx = 1:length(file_list)
    [~, name, ~] = fileparts(file_list(idx).name);
    filename = fullfile(config.input_folder, file_list(idx).name);
    
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('处理图像 %d/%d: %s\n', idx, length(file_list), name);
    fprintf('%s\n', repmat('=', 1, 60));
    fprintf('当前参数: 分块=%d×%d, 阶数=%d, eta=%.3f, 迭代=%d\n', ...
        legendre_param.num_blocks_m, legendre_param.num_blocks_n, ...
        legendre_param.poly_degree, energy_param.eta, energy_param.max_iter);
    
    %% 加载图像
    img = imread(filename);
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    g = im2double(img);
    
    [M, N] = size(g);
    fprintf('图像尺寸: %d × %d\n', M, N);
    
    %% 热辐射偏置场拟合
    tic;
    
    if use_energy_minimization
        [fb, b, cost_history] = minimize_energy_legendre(g, ...
            'eta', energy_param.eta, ...
            'max_iter', energy_param.max_iter, ...
            'tol', energy_param.tol, ...
            'display', energy_param.display, ...
            'poly_degree', legendre_param.poly_degree, ...
            'num_blocks_m', legendre_param.num_blocks_m, ...
            'num_blocks_n', legendre_param.num_blocks_n);
    else
        % b = legendre_surface_fitting_adaptive(g, legendre_param);
        % fb = g - b;
        % cost_history = [];
        % --- 获取原始偏置场 ---
    b = legendre_surface_fitting_adaptive(g, legendre_param);
    
    % --- 关键优化：对偏置场 b 进行轻微平滑 (就地修改) ---
    % 定义一个非常轻微的高斯平滑核
    sigma = 17; % 平滑强度，必须很小！从0.3, 0.5开始尝试
    kernel_size = 17;% 核大小，奇数，如3, 5
    
    % 创建高斯核并应用平滑
    gaussian_kernel = fspecial('gaussian', kernel_size, sigma);
    b = imfilter(b, gaussian_kernel, 'conv', 'replicate'); % 注意：这里直接覆盖了 b
    
    % --- 计算校正后图像 (fb) ---
    fb = g - b; % 使用平滑后的 b
    
    elapsed_time = toc;
    cost_history = [];
    end
    
    elapsed_time = toc;
    fprintf('处理耗时: %.3f 秒\n', elapsed_time);
    
    %% 后处理：归一化
    fb = norm01(fb);
    b = norm01(b);
    
    %% 计算图像质量指标
    % 参考图像（如果存在）
    ref_filename = fullfile(config.clearer_folder, [name, '.png']);
    if exist(ref_filename, 'file')
        ref_img = imread(ref_filename);
        if size(ref_img, 3) == 3
            ref_img = rgb2gray(ref_img);
        end
        ref_img = im2double(ref_img);
        ref_img = norm01(ref_img);
        psnr_ref = calculatePSNR(ref_img, fb);
        ssim_ref = calculateSSIM(ref_img, fb);
    else
        psnr_ref = NaN;
        ssim_ref = NaN;
    end
    
    % 无参考指标
    psnr_val = calculatePSNR(g, fb);
    ssim_val = calculateSSIM(g, fb);
    eog_val = EOG(fb);
    niqe_val = niqe(fb);
    nisa_val = final(fb);
    niqe_orig = niqe(g);
    nisa_orig = final(g);
    
    %% 显示结果
    fprintf('\n--- 质量指标 ---\n');
    fprintf('PSNR: %.4f dB\n', psnr_val);
    fprintf('SSIM: %.4f\n', ssim_val);
    fprintf('EOG:  %.4f\n', eog_val);
    fprintf('NIQE: %.4f (原始: %.4f)\n', niqe_val, niqe_orig);
    fprintf('NISA: %.4f (原始: %.4f)\n', nisa_val, nisa_orig);
    
    %% 保存结果
    if save_results
        % 保存校正图像
        imwrite(fb, fullfile(config.output_folder, [name, '_corrected.png']));
        % 保存偏置场
        imwrite(b, fullfile(config.output_folder, [name, '_bias.png']));
        
        % 保存3D曲面图
        fig = figure('Visible', 'off');
        surf(b);
        shading interp;
        colormap jet;
        title(sprintf('热辐射偏置场 3D (%s)', name));
        saveas(fig, fullfile(config.output_folder, [name, '_bias_3D.png']));
        close(fig);
        
        fprintf('结果已保存到: %s\n', config.output_folder);
    end
    
    %% 显示结果图
    if show_results
        figure('Name', sprintf('处理结果: %s', name), 'Position', [100, 100, 1400, 500]);
        
        subplot(1,4,1);
        imshow(g);
        title('原始图像');
        xlabel(sprintf('NIQE=%.2f, NISA=%.2f', niqe_orig, nisa_orig));
        
        subplot(1,4,2);
        imshow(fb);
        title('校正后图像');
        xlabel(sprintf('PSNR=%.2f, SSIM=%.3f', psnr_val, ssim_val));
        
        subplot(1,4,3);
        imshow(b);
        title('热辐射偏置场');
        
        subplot(1,4,4);
        imshow(g - fb, []);
        title('残差图像');
        colormap(gca, 'jet');
        
        sgtitle(sprintf('%s (勒让德阶数=%d, 耗时=%.2fs)', name, legendre_param.poly_degree, elapsed_time));
    end
    
    %% 存储结果
    results{idx, 1} = name;
    results{idx, 2} = psnr_val;
    results{idx, 3} = ssim_val;
    results{idx, 4} = eog_val;
    results{idx, 5} = niqe_val;
    results{idx, 6} = nisa_val;
    results{idx, 7} = niqe_orig;
    results{idx, 8} = nisa_orig;
    results{idx, 9} = psnr_ref;
    results{idx, 10} = ssim_ref;
    results{idx, 11} = elapsed_time;
    results{idx, 12} = legendre_param.poly_degree;
end

%% 保存结果表格
results_table = cell2table(results, 'VariableNames', result_names);
writetable(results_table, fullfile(config.output_folder, 'results_summary.csv'));
fprintf('\n结果汇总已保存到: %s\n', fullfile(config.output_folder, 'results_summary.csv'));

%% 显示汇总
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('处理完成！结果汇总:\n');
fprintf('%s\n', repmat('=', 1, 80));
fprintf('平均 PSNR: %.4f dB\n', mean([results{:, 2}]));
fprintf('平均 SSIM: %.4f\n', mean([results{:, 3}]));
fprintf('平均 EOG:  %.4f\n', mean([results{:, 4}]));
fprintf('平均 NIQE: %.4f\n', mean([results{:, 5}]));
fprintf('平均 NISA: %.4f\n', mean([results{:, 6}]));
fprintf('平均耗时: %.3f 秒\n', mean([results{:, 11}]));
disp(results_table);
fprintf('程序运行完成！\n');