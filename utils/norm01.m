function M = norm01(Mat)

M = (Mat-min(Mat(:)))./(max(Mat(:))-min(Mat(:)));

end