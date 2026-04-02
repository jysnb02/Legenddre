function  saveas_hard(I,dir)
if ~isfolder(fileparts(dir))
    mkdir(fileparts(dir));
end
    saveas(I,dir);
end

