function batchFormatDb(filepath)

nf=length(filepath);

for n = 1:nf
    commandtext=['formatdb.exe -i ' filepath{n,1} ' -o T -p F'];
    system(commandtext);
end

