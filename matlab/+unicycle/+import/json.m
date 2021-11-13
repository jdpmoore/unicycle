function obj = json(fname)
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
obj = jsondecode(str);
end