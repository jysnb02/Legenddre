function results = process_single_image(image_path, out_folder, varargin)
% results = process_single_image(image_path, out_folder, ...)
% 简单封装，处理单张图像并保存结果。
% 可选参数 (name,value):
%   'poly_degree' (default 3)
%   'num_blocks_m' (default 2)
%   'num_blocks_n' (default 2)
%   'use_energy_minimization' (default true)
%   'eta' (energy minimization eta, default 0.01)
%   'max_iter' (default 5)
%   'tol' (default 1e-6)
%   'save_results' (default true)
%   'show_results' (default false)

% parse inputs
p = inputParser;
addRequired(p,'image_path',@ischar);
addRequired(p,'out_folder',@ischar);
addParameter(p,'poly_degree',3,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'num_blocks_m',2,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'num_blocks_n',2,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'use_energy_minimization',true,@islogical);
addParameter(p,'eta',0.01,@isnumeric);
addParameter(p,'max_iter',5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'tol',1e-6,@isnumeric);
addParameter(p,'save_results',true,@islogical);
addParameter(p,'show_results',false,@islogical);
parse(p,image_path,out_folder,varargin{:});
opts = p.Results;

% prepare
if ~exist(opts.out_folder,'dir')
    mkdir(opts.out_folder);
end

% read image
img = imread(opts.image_path);
if size(img,3) == 3
    img = rgb2gray(img);
end
g = im2double(img);

% run processing
tic;
if opts.use_energy_minimization
    [fb, b, cost_history] = minimize_energy_legendre(g, ...
        'eta', opts.eta, ...
        'max_iter', opts.max_iter, ...
        'tol', opts.tol, ...
        'display', false, ...
        'poly_degree', opts.poly_degree, ...
        'num_blocks_m', opts.num_blocks_m, ...
        'num_blocks_n', opts.num_blocks_n);
else
    param.num_blocks_m = opts.num_blocks_m;
    param.num_blocks_n = opts.num_blocks_n;
    param.poly_degree = opts.poly_degree;
    b = legendre_surface_fitting_adaptive(g, param);
    fb = g - b;
    cost_history = [];
end
elapsed = toc;

% normalize (same as main_legendre)
fbn = norm01(fb);
bn = norm01(b);

% metrics (if functions available)
try
    psnr_val = calculatePSNR(g, fbn);
catch
    psnr_val = NaN;
end
try
    ssim_val = calculateSSIM(g, fbn);
catch
    ssim_val = NaN;
end
try
    eog_val = EOG(fbn);
catch
    eog_val = NaN;
end

% save outputs
[~, name, ~] = fileparts(opts.image_path);
if opts.save_results
    imwrite(fbn, fullfile(opts.out_folder, [name, '_corrected.png']));
imwrite(bn, fullfile(opts.out_folder, [name, '_bias.png']));
    % save 3D bias figure
    fig = figure('Visible','off'); surf(bn); shading interp; colormap jet;
    saveas(fig, fullfile(opts.out_folder, [name, '_bias_3D.png']));
    close(fig);
    % save comparison
    fig = figure('Visible','off');
    subplot(1,3,1); imshow(g); title('原始');
    subplot(1,3,2); imshow(fbn); title('校正后');
    subplot(1,3,3); imshow(bn); title('偏置场');
    saveas(fig, fullfile(opts.out_folder, [name, '_comparison.png']));
    close(fig);
end

% package results
results.original = g;
results.corrected = fbn;
results.bias = bn;
results.psnr = psnr_val;
results.ssim = ssim_val;
results.eog = eog_val;
results.cost_history = cost_history;
results.elapsed = elapsed;
fprintf('Processed %s in %.3f s. PSNR=%.4f, SSIM=%.4f\n', name, elapsed, psnr_val, ssim_val);
end