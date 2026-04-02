function imwrite_hard(I,path)
if ~isfolder(fileparts(path))
    mkdir(fileparts(path));
end
imwrite(I,path);
end

