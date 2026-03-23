%% ============================================================
% 勒让德多项式热辐射偏置场校正主程序
% 功能：读取图像，进行勒让德曲面拟合，去除热辐射偏置场
% 支持自适应阶数选择（1-5阶）
% ============================================================

clc;
clear;
close all;
addpath('./utils/');

%% 1. 用户设置 - 请在这里修改参数
% ============================================================
% 输入图像路径（修改为你的图像路径）
input_folder = 'E:\test_legenddre\仿真图\原图';
output_folder = 'E:\test_legenddre\output_legendre\';

% 分块参数
param.num_blocks_m = 2;      % 垂直方向分块数
param.num_blocks_n = 2;      % 水平方向分块数
param.poly_degree = 3;       % 多项式阶数 (1,2,3,4,5)
param.save_blocks = true;    % 是否保存中间结果
param.save_path = output_folder;

% 能量最小化参数
energy_params.eta = 0.01;     % 梯度保真项权重
energy_params.max_iter = 5;   % 最大迭代次数
energy_params.tol = 1e-6;     % 收敛容差
energy_params.display = true; % 显示迭代信息

% 是否使用能量最小化（true: 交替优化, false: 直接拟合）
use_energy_minimization = true;

% 是否保存结果图像
save_results = true;

% 是否显示结果图
show_results = true;

%% 2. 创建输出文件夹
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% 3. 获取输入图像列表
file_list = dir(fullfile(input_folder, '*.png'));
if isempty(file_list)
    file_list = dir(fullfile(input_folder, '*.jpg'));
end
if isempty(file_list)
    file_list = dir(fullfile(input_folder, '*.bmp'));
end
if isempty(file_list)
    file_list = dir(fullfile(input_folder, '*.tif'));
end

if isempty(file_list)
    error('未找到图像文件，请检查路径: %s', input_folder);
end

fprintf('找到 %d 张图像\n\n', length(file_list));

%% 4. 初始化结果存储
results = cell(length(file_list), 12);
result_names = {'图像名', 'PSNR', 'SSIM', 'EOG', 'NIQE', 'NISA', ...
                'NIQE_orig', 'NISA_orig', 'PSNR_ref', 'SSIM_ref', ...
                '耗时(s)', '阶数'};

%% 5. 主循环：处理每一张图像
for idx = 1:length(file_list)
    [~, name, ~] = fileparts(file_list(idx).name);
    filename = fullfile(input_folder, file_list(idx).name);
    
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('处理图像 %d/%d: %s\n', idx, length(file_list), name);
    fprintf('%s\n', repmat('=', 1, 60));
    
    %% 5.1 加载图像
    img = imread(filename);
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    g = im2double(img);
    
    % 记录原始图像尺寸
    [M, N] = size(g);
    fprintf('图像尺寸: %d × %d\n', M, N);
    
    %% 5.2 热辐射偏置场拟合
    tic;
    
    if use_energy_minimization
        % 使用能量最小化（交替优化）
        [fb, b, cost_history] = minimize_energy_legendre(g, ...
            'eta', energy_params.eta, ...
            'max_iter', energy_params.max_iter, ...
            'tol', energy_params.tol, ...
            'display', energy_params.display, ...
            'poly_degree', param.poly_degree, ...
            'num_blocks_m', param.num_blocks_m, ...
            'num_blocks_n', param.num_blocks_n);
    else
        % 直接使用勒让德曲面拟合
        b = legendre_surface_fitting_adaptive(g, param);
        fb = g - b;
        cost_history = [];
    end
    
    elapsed_time = toc;
    fprintf('处理耗时: %.3f 秒\n', elapsed_time);
    
    %% 5.3 后处理：归一化到 [0,1]
    fb = norm01(fb);
    b = norm01(b);
    
    %% 5.4 计算图像质量指标（如果有参考图像）
    % 尝试加载对应的清晰图像（假设文件名相同，放在clearer文件夹）
    clearer_folder = fullfile(input_folder, '清晰图');
    ref_filename = fullfile(clearer_folder, [name, '.png']);
    
    if exist(ref_filename, 'file')
        ref_img = imread(ref_filename);
        if size(ref_img, 3) == 3
            ref_img = rgb2gray(ref_img);
        end
        ref_img = im2double(ref_img);
        ref_img = norm01(ref_img);
        
        % 计算参考指标
        psnr_ref = calculatePSNR(ref_img, fb);
        ssim_ref = calculateSSIM(ref_img, fb);
        fprintf('参考图像 PSNR: %.4f dB, SSIM: %.4f\n', psnr_ref, ssim_ref);
    else
        psnr_ref = NaN;
        ssim_ref = NaN;
        fprintf('未找到参考图像: %s\n', ref_filename);
    end
    
    % 计算无参考指标
    psnr_val = calculatePSNR(g, fb);
    ssim_val = calculateSSIM(g, fb);
    eog_val = EOG(fb);
    niqe_val = niqe(fb);
    nisa_val = final(fb);
    
    % 原始图像指标
    niqe_orig = niqe(g);
    nisa_orig = final(g);
    
    %% 5.5 显示结果
    fprintf('\n--- 质量指标 ---\n');
    fprintf('PSNR: %.4f dB\n', psnr_val);
    fprintf('SSIM: %.4f\n', ssim_val);
    fprintf('EOG:  %.4f\n', eog_val);
    fprintf('NIQE: %.4f (原始: %.4f)\n', niqe_val, niqe_orig);
    fprintf('NISA: %.4f (原始: %.4f)\n', nisa_val, nisa_orig);
    
    %% 5.6 保存结果
    if save_results
        % 保存校正后的图像
        imwrite(fb, fullfile(output_folder, [name, '_corrected.png']));
        
        % 保存偏置场图像
        imwrite(b, fullfile(output_folder, [name, '_bias.png']));
        
        % 保存3D曲面图
        fig = figure('Visible', 'off');
        surf(b);
        shading interp;
        colormap jet;
        title(sprintf('热辐射偏置场 3D (%s)', name));
        saveas(fig, fullfile(output_folder, [name, '_bias_3D.png']));
        close(fig);
        
        % 保存对比图
        fig = figure('Visible', 'off');
        subplot(2,3,1); imshow(g); title('原始图像');
        subplot(2,3,2); imshow(fb); title('校正后图像');
        subplot(2,3,3); imshow(b); title('偏置场');
        subplot(2,3,4); imshow(g - fb, []); title('残差');
        subplot(2,3,5); imshow(fb + b); title('重建图像');
        if ~isempty(cost_history)
            subplot(2,3,6); semilogy(cost_history); title('收敛曲线'); xlabel('迭代'); ylabel('成本');
        end
        sgtitle(sprintf('%s (阶数=%d, PSNR=%.2f, SSIM=%.3f)', name, param.poly_degree, psnr_val, ssim_val));
        saveas(fig, fullfile(output_folder, [name, '_comparison.png']));
        close(fig);
        
        fprintf('结果已保存到: %s\n', output_folder);
    end
    
    %% 5.7 显示结果图
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
        
        sgtitle(sprintf('%s (勒让德阶数=%d, 耗时=%.2fs)', name, param.poly_degree, elapsed_time));
    end
    
    %% 5.8 存储结果
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
    results{idx, 12} = param.poly_degree;
end

%% 6. 保存结果表格
% ============================================================
results_table = cell2table(results, 'VariableNames', result_names);
writetable(results_table, fullfile(output_folder, 'results_summary.csv'));
fprintf('\n结果汇总已保存到: %s\n', fullfile(output_folder, 'results_summary.csv'));

%% 7. 显示结果汇总
% ============================================================
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('处理完成！结果汇总:\n');
fprintf('%s\n', repmat('=', 1, 80));

% 计算平均指标
avg_psnr = mean([results{:, 2}]);
avg_ssim = mean([results{:, 3}]);
avg_eog = mean([results{:, 4}]);
avg_niqe = mean([results{:, 5}]);
avg_nisa = mean([results{:, 6}]);
avg_time = mean([results{:, 11}]);

fprintf('平均 PSNR:  %.4f dB\n', avg_psnr);
fprintf('平均 SSIM:  %.4f\n', avg_ssim);
fprintf('平均 EOG:   %.4f\n', avg_eog);
fprintf('平均 NIQE:  %.4f\n', avg_niqe);
fprintf('平均 NISA:  %.4f\n', avg_nisa);
fprintf('平均耗时:   %.3f 秒\n', avg_time);

% 显示结果表格
disp(results_table);

fprintf('\n所有结果已保存到: %s\n', output_folder);
fprintf('程序运行完成！\n');