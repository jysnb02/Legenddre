function M = norm_11(Mat)

M = 2.*(Mat-min(Mat(:)))./(max(Mat(:))-min(Mat(:)))-1;

end