function seq = fasTrim(seq_origin, seq_delete)

seq = seq_origin;

totalnum = [];
for n = 1:length(seq_delete)
    pos1 = strfind(seq_delete{n,1}, '_');
    pos2 = strfind(seq_delete{n,1}, '=');
    temp = char(seq_delete{n,1});
    if ~isempty(pos2)
        num = str2double(temp(pos2(end)+1:end));
        num = uint32(num);
    else
        num = str2double(temp(pos1+1:end));
        num = uint32(num);
    end
    totalnum = [totalnum 2*num-1 2*num];
end

seq(totalnum) = [];

