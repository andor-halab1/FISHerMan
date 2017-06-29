function [Header,Sequence]...
    =abundantrnaParse(cdnaHeader,cdnaSequence,ncrnaHeader,ncrnaSequence,seqData,params)

% params = struct('species','Mouse','verbose',1,...
%     'percent',0.001,...
%     'keys',{':rRNA',':Mt_rRNA',':tRNA',':Mt_tRNA'});

if params(1).verbose
    disp('generating abundant rna database files for Blast');
end

% pick highly abundant rna based on seqData
transcriptID = {};
if ~isempty(seqData)
    seqDataTable = cell2table(seqData);
    seqDataTable = sortrows(seqDataTable,[3 1],{'descend','ascend'});
    seqData = table2cell(seqDataTable);
    geneNum = floor(size(seqData,1)*params(1).percent);
    transcriptID(end+1:end+geneNum,1) = seqData(1:geneNum,1);
end

% make sure rRNA and tRNA are included
for i = 1:length(params)
    pos{i,1} = regexp(ncrnaHeader, params(i).keys, 'end');
end

for n = 1:length(ncrnaHeader)
    if ~(isempty(pos{1,1}{n,1}) && isempty(pos{2,1}{n,1})...
            && isempty(pos{3,1}{n,1}) && isempty(pos{4,1}{n,1}))
        temp=regexp(ncrnaHeader{n,1}, ':');
        transcriptID{end+1,1} = ncrnaHeader{n,1}(1:temp(1)-1);
    end
end
transcriptID = unique(transcriptID);

Header = vertcat(cdnaHeader,ncrnaHeader);
Sequence = vertcat(cdnaSequence,ncrnaSequence);
[Header, Sequence] = pickExpressedSeq(transcriptID, Header, Sequence);

abundantrna = [params(1).species '.abundantrna.fas'];
if exist(abundantrna, 'file')
    delete(abundantrna);
end
fastawrite(abundantrna, Header, Sequence);

abundantrnaDb = [params(1).species '.abundantrnaDb.fas'];
% MatLab's use of blastlocal requires short entry name
simpleHeader = Header;
for n = 1:length(Header)
    pos = regexp(Header{n,1}, ':');
    simpleHeader{n,1} = Header{n,1}(1:pos(1)-1);
end
if exist(abundantrnaDb, 'file')
    delete([abundantrnaDb '*']);
end
fastawrite(abundantrnaDb, simpleHeader, Sequence);
blastformat('Inputdb', abundantrnaDb,...
    'FormatArgs', '-o T -p F');
