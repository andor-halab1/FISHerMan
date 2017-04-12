function seq = fasRead(filepath, len)

if nargin == 1
    len = 480000;
end

fid = fopen(filepath,'r');
seq = textscan(fid, '%s', len, 'Delimiter', '\n');
seq = seq{1,1};
fclose(fid);

end

