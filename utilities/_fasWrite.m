function fasWrite(filepath, seq)

fid = fopen(filepath,'a');
for n = 2:2:length(seq)
    fprintf(fid, '%s\n', strcat(seq{n-1,1}, '=', num2str(n/2)));
    fprintf(fid, '%s\n', seq{n,1});
end
fclose(fid);

end