function [Header,Sequence]=ncrnaParse(ncrna,seqData,trna,params)

% params = struct('species','Mouse','verbose',1,...
%     'dir1','C:\FISHerMan\Db\Mouse38.ncrna.fa',...
%     'tRNA',1,'dirT','C:\FISHerMan\Db\Mouse.trna.fas',...
%     'keys',{'ncrna','gene:\S*','gene_biotype:\S*'});

if params(1).verbose
    disp('reading the ncrna data file');
end

[Header, Sequence] = fastaread(ncrna);
Header = Header';
Sequence = Sequence';

if params(1).verbose
    disp('  trimming fasta headers');
end

for i = 1:length(params)
    [pos1{i,1}, pos2{i,1}] = regexp(Header, params(i).keys, 'start', 'end');
end

for n = 1:length(Header)
    temp = Header{n,1};
    temp1 = temp(1:pos1{1,1}{n,1}-2);
    temp2 = temp(pos1{2,1}{n,1}+5:pos2{2,1}{n,1});
    temp3 = temp(pos1{3,1}{n,1}+13:pos2{3,1}{n,1});
    
    if isempty(temp1)
        disp('missing transcript ID');
    elseif strfind(temp1,'.')
        temp1pos=strfind(temp1,'.');
        temp1=temp1(1:temp1pos(1)-1);
    end
    if isempty(temp2)
        disp('missing gene ID');
    elseif strfind(temp2,'.')
        temp2pos=strfind(temp2,'.');
        temp2=temp2(1:temp2pos(1)-1);
    end
    if isempty(temp3)
        disp('missing gene type');
    elseif strfind(temp3,'.')
        temp3pos=strfind(temp3,'.');
        temp3=temp3(1:temp3pos(1)-1);
    end
    
    Header{n,1} = strcat(temp1, ':', temp2, ':', temp3);
end

if ~isempty(seqData)
    if params(1).verbose
        disp('  picking expressed sequences according to RNA-seq data');
    end
    [Header, Sequence] = pickExpressedSeq(seqData, Header, Sequence);
end

if ~isempty(trna)
    if params(1).verbose
        disp('  appending additional tRNA sequences to ncrna data file');
    end
    [trnaHeader, trnaSequence] = fastaread(trna);
    trnaHeader = trnaHeader';
    trnaSequence = trnaSequence';
    for n = 1:length(trnaHeader)
        Header{end+1,1}=['ENSRNAT' num2str(n) ':ENSRNAG' num2str(n) ':tRNA'];
        Sequence{end+1,1}=trnaSequence{n,1};
    end
end

if params(1).verbose
    disp('saving the ncrna fasta file');
end

ncrna = [params(1).species '.ncrna.fas'];
if exist(ncrna, 'file')
    delete(ncrna);
end
fastawrite(ncrna, Header, Sequence);

if params(1).verbose
    disp('generating ncrna database files for Blast');
end

ncrnaDb = [params(1).species '.ncrnaDb.fas'];
% MatLab's use of blastlocal requires short entry names
simpleHeader = Header;
for n = 1:length(Header)
    pos = regexp(Header{n,1}, ':');
    simpleHeader{n,1} = Header{n,1}(1:pos(1)-1);
end
if exist(ncrnaDb, 'file')
    delete([ncrnaDb '*']);
end
fastawrite(ncrnaDb, simpleHeader, Sequence);
blastformat('Inputdb', ncrnaDb,...
    'FormatArgs', '-o T -p F');

