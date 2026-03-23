function score = EOG(Image)
% 转换为灰度图
if size(Image, 3) == 3
    Image = rgb2gray(Image);
end
Image = double(Image);
[rows, cols] = size(Image);

% 计算梯度（使用简单差分）
gradient_x = diff(Image, 1, 2);
gradient_x = [gradient_x, zeros(rows, 1)]; % 补边
gradient_y = diff(Image, 1, 1);
gradient_y = [gradient_y; zeros(1, cols)];

% 计算梯度能量
gradient_energy = gradient_x.^2 + gradient_y.^2;
score = sum(gradient_energy(:));
end