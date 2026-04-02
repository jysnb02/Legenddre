function psnr_value = calculatePSNR(original, distorted)
    % 将图像转换为double类型进行计算
    original = double(original);
    distorted = double(distorted);
    
    % 计算均方误差 (MSE)
    mse = mean((original(:) - distorted(:)).^2);
    
    % 设定最大像素值（对于8位图像为255）
    max_pixel = 255.0;
    
    % 计算PSNR
    if mse == 0
        psnr_value = 100; % 如果两图完全相同，MSE为0，PSNR为无穷大，这里设为一个较大值
    else
        psnr_value = 10 * log10((max_pixel^2) / mse);
    end
end