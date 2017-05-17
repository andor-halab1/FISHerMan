function params=readParameters(species,parameters)

% species = 'Mouse';
% parameters = 'Mouse.parameters.xml';

params = struct('rnaSeq','','cdna','','ncrna','',...
    'abundantrna','','transcriptList','','OligoArray','',...
    'oligos','');

tree = xmlread(parameters);

try
    parameters = tree.getElementsByTagName('parameters');
    parameters = parameters.item(0);
    
    general = parameters.getElementsByTagName('general');
    general = general.item(0);
    
    rnaSeq = parameters.getElementsByTagName('rnaSeq');
    rnaSeq = rnaSeq.item(0);
    
    cdna = parameters.getElementsByTagName('cdna');
    cdna = cdna.item(0);
    
    ncrna = parameters.getElementsByTagName('ncrna');
    ncrna = ncrna.item(0);
    
    abundantrna = parameters.getElementsByTagName('abundantrna');
    abundantrna = abundantrna.item(0);
    
    transcriptList = parameters.getElementsByTagName('transcriptList');
    transcriptList = transcriptList.item(0);
    
    OligoArray = parameters.getElementsByTagName('OligoArray');
    OligoArray = OligoArray.item(0);
    
    oligos = parameters.getElementsByTagName('oligos');
    oligos = oligos.item(0);
catch
    disp('error occurred at level 1');
end

%%
try
    verbose = general.getElementsByTagName('verbose');
    verbose = verbose.item(0);
catch
    disp('error occurred at general');
end

%%
try
    key1 = rnaSeq.getElementsByTagName('key1');
    key1 = key1.item(0);
    
    key2 = rnaSeq.getElementsByTagName('key2');
    key2 = key2.item(0);
    
    thres = rnaSeq.getElementsByTagName('thres');
    thres = thres.item(0);
catch
    disp('error occurred at rnaSeq');
end

params.rnaSeq = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'keys',{char(key1.getFirstChild.getData),char(key2.getFirstChild.getData)},...
    'thres',str2double(thres.getFirstChild.getData));

%%
try
    key1 = cdna.getElementsByTagName('key1');
    key1 = key1.item(0);
    
    key2 = cdna.getElementsByTagName('key2');
    key2 = key2.item(0);
catch
    disp('error occurred at cdna');
end

params.cdna = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'keys',{char(key1.getFirstChild.getData),char(key2.getFirstChild.getData)});

%%
try
    key1 = ncrna.getElementsByTagName('key1');
    key1 = key1.item(0);
    
    key2 = ncrna.getElementsByTagName('key2');
    key2 = key2.item(0);
    
    key3 = ncrna.getElementsByTagName('key3');
    key3 = key3.item(0);
catch
    disp('error occurred at ncrna');
end

params.ncrna = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'keys',{char(key1.getFirstChild.getData),...
    char(key2.getFirstChild.getData),char(key3.getFirstChild.getData)});

%%
try
    key1 = abundantrna.getElementsByTagName('key1');
    key1 = key1.item(0);
    
    key2 = abundantrna.getElementsByTagName('key2');
    key2 = key2.item(0);
    
    key3 = abundantrna.getElementsByTagName('key3');
    key3 = key3.item(0);
    
    key4 = abundantrna.getElementsByTagName('key4');
    key4 = key4.item(0);
    
    key5 = abundantrna.getElementsByTagName('key5');
    key5 = key5.item(0);
catch
    disp('error occurred at abundantrna');
end

params.abundantrna = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'keys',{char(key1.getFirstChild.getData),...
    char(key2.getFirstChild.getData),char(key3.getFirstChild.getData),...
    char(key4.getFirstChild.getData),char(key5.getFirstChild.getData)});

%%
try
    len = transcriptList.getElementsByTagName('length');
    len = len.item(0);
    
    num = transcriptList.getElementsByTagName('number');
    num = num.item(0);
catch
    disp('error occurred at transcriptList');
end

params.transcriptList = struct('species',species,...
    'verbose',str2double(verbose.getFirstChild.getData),...
    'length',str2double(len.getFirstChild.getData),...
    'number',str2double(num.getFirstChild.getData));

%%

%%
try
    key1 = oligos.getElementsByTagName('key1');
    key1 = key1.item(0);
    
    key2 = oligos.getElementsByTagName('key2');
    key2 = key2.item(0);
    
    num = oligos.getElementsByTagName('number');
    num = num.item(0);
    
    thres = oligos.getElementsByTagName('thres');
    thres = thres.item(0);
    
    querySize = oligos.getElementsByTagName('querySize');
    querySize = querySize.item(0);
    
    DbSize = oligos.getElementsByTagName('DbSize');
    DbSize = DbSize.item(0);
    
    seqNum = oligos.getElementsByTagName('seqNum');
    seqNum = seqNum.item(0);
    
    blastArgs = oligos.getElementsByTagName('blastArgs');
    blastArgs = blastArgs.item(0);
    
    parallel = oligos.getElementsByTagName('parallel');
    parallel = parallel.item(0);
catch
    disp('error occurred at oligos');
end

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

