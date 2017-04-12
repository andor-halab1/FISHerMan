seq_origin = fasRead('C:\Users\frank\Documents\MATLAB\Probe_Gen\probes_egp1c_e.fas');
seq_origin = seq_origin(1:2:end);

totalnum = [];
for n = 1:length(seq_origin)
    pos1 = strfind(seq_origin{n,1}, '_');
    pos2 = strfind(seq_origin{n,1}, '=');
    temp = char(seq_origin{n,1});
    num = str2double(temp(pos1+1:pos2(1)-1));
    num = uint32(num);
    totalnum = [totalnum num];
end
totalnum = totalnum';
totalnum = sort(totalnum);

totalnum = ceil(double(totalnum)/96.0);
totalnum = uint32(totalnum);
totalnum_unique = unique(totalnum, 'stable');

name_occurrence = [];
for n = 1:length(totalnum_unique)
    temp = find(totalnum == totalnum_unique(n));
    temp = length(temp);
    name_occurrence = [name_occurrence temp];
end
name_occurrence = name_occurrence';

