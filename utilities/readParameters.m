function params=readParameters(species,parameters)

% species = 'Mouse';
% parameters = 'Mouse.parameters.xml';

params = struct('rnaSeq','','cdna','','ncrna','',...
    'abundantrna','','transcriptList','','OligoArray','',...
    'oligos','');

tree = xmlread(parameters);

parameters=xmlParse(tree, 'tree', 'parameters');
general=xmlParse(parameters, 'parameters', 'general');
rnaSeq=xmlParse(parameters, 'parameters', 'rnaSeq');
cdna=xmlParse(parameters, 'parameters', 'cdna');
ncrna=xmlParse(parameters, 'parameters', 'ncrna');
abundantrna=xmlParse(parameters, 'parameters', 'abundantrna');
transcriptList=xmlParse(parameters, 'parameters', 'transcriptList');
OligoArray=xmlParse(parameters, 'parameters', 'OligoArray');
oligos=xmlParse(parameters, 'parameters', 'oligos');

%%
verbose=xmlParse(general, 'general', 'verbose');

%%
key1=xmlParse(rnaSeq, 'rnaSeq', 'key1');
key2=xmlParse(rnaSeq, 'rnaSeq', 'key2');
thres=xmlParse(rnaSeq, 'rnaSeq', 'thres');

params.rnaSeq = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'keys',{char(key1.getFirstChild.getData),char(key2.getFirstChild.getData)},...
    'thres',str2double(thres.getFirstChild.getData));

%%
key1=xmlParse(cdna, 'cdna', 'key1');
key2=xmlParse(cdna, 'cdna', 'key2');

params.cdna = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'keys',{char(key1.getFirstChild.getData),char(key2.getFirstChild.getData)});

%%
key1=xmlParse(ncrna, 'ncrna', 'key1');
key2=xmlParse(ncrna, 'ncrna', 'key2');
key3=xmlParse(ncrna, 'ncrna', 'key3');

params.ncrna = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'keys',{char(key1.getFirstChild.getData),...
    char(key2.getFirstChild.getData),char(key3.getFirstChild.getData)});

%%
percent=xmlParse(abundantrna, 'abundantrna', 'percent');
key1=xmlParse(abundantrna, 'abundantrna', 'key1');
key2=xmlParse(abundantrna, 'abundantrna', 'key2');
key3=xmlParse(abundantrna, 'abundantrna', 'key3');
key4=xmlParse(abundantrna, 'abundantrna', 'key4');
key5=xmlParse(abundantrna, 'abundantrna', 'key5');

params.abundantrna = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'percent',str2double(percent.getFirstChild.getData),...
    'keys',{char(key1.getFirstChild.getData),...
    char(key2.getFirstChild.getData),char(key3.getFirstChild.getData),...
    char(key4.getFirstChild.getData),char(key5.getFirstChild.getData)});

%%
len=xmlParse(transcriptList, 'transcriptList', 'length');
num=xmlParse(transcriptList, 'transcriptList', 'number');

params.transcriptList = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'length',str2double(len.getFirstChild.getData),...
    'number',str2double(num.getFirstChild.getData));

%%

%%
key1=xmlParse(oligos, 'oligos', 'key1');
key2=xmlParse(oligos, 'oligos', 'key2');
num=xmlParse(oligos, 'oligos', 'number');
thres=xmlParse(oligos, 'oligos', 'thres');
querySize=xmlParse(oligos, 'oligos', 'querySize');
DbSize=xmlParse(oligos, 'oligos', 'DbSize');
seqNum=xmlParse(oligos, 'oligos', 'seqNum');
blastArgs=xmlParse(oligos, 'oligos', 'blastArgs');
parallel=xmlParse(oligos, 'oligos', 'parallel');

params.oligos = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'keys',{char(key1.getFirstChild.getData),char(key2.getFirstChild.getData)},...
    'number',str2double(num.getFirstChild.getData),...
    'thres',str2double(thres.getFirstChild.getData),...
    'querySize',str2double(querySize.getFirstChild.getData),...
    'DbSize',str2double(DbSize.getFirstChild.getData),...
    'seqNum',str2double(seqNum.getFirstChild.getData),...
    'blastArgs',char(blastArgs.getFirstChild.getData),...
    'parallel',str2double(parallel.getFirstChild.getData));

