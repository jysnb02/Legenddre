function ssim_value = calculateSSIM(img1, img2)
    % 将图像转换为double类型
    img1 = double(img1);
    img2 = double(img2);
    
    % 常数定义，防止分母为零
    K1 = 0.01;
    K2 = 0.03;
    L = 255; % 像素值的动态范围（8位图像为255）
    C1 = (K1 * L)^2;
    C2 = (K2 * L)^2;
    
    % 使用一个小窗口（例如11x11的高斯加权函数）来计算局部统计量
    % 这里为简化，使用整个图像计算全局SSIM
    mu1 = mean(img1(:));
    mu2 = mean(img2(:));
    sigma1 = std(img1(:), 1); % 使用总体标准差
    sigma2 = std(img2(:), 1);
    sigma12 = mean((img1(:) - mu1) .* (img2(:) - mu2));
    
    % 计算SSIM
    numerator = (2 * mu1 * mu2 + C1) * (2 * sigma12 + C2);
    denominator = (mu1^2 + mu2^2 + C1) * (sigma1^2 + sigma2^2 + C2);
    ssim_value = numerator / denominator;
end