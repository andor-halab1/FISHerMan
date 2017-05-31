function [Header,Sequence]=cdnaParse(cdna,seqData,params)

% params = struct('species','Mouse','verbose',1,...
%     'dir1','C:\FISHerMan\Db\Mouse38.cdna.fa',...
%     'keys',{'cdna','gene:\S*'});

if params(1).verbose
    disp('reading the cdna data file');
end

[Header, Sequence] = fastaread(cdna);
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
    
    Header{n,1} = strcat(temp1, ':', temp2);
end

if ~isempty(seqData)
    if params(1).verbose
        disp('  picking expressed sequences according to RNA-seq data');
    end
    [Header, Sequence] = pickExpressedSeq(seqData, Header, Sequence);
end

if params(1).verbose
    disp('saving the cdna fasta file');
end

cdna = [params(1).species '.cdna.fas'];
if exist(cdna, 'file')
    delete(cdna);
end
fastawrite(cdna, Header, Sequence);

if params(1).verbose
    disp('generating cdna database files for Blast');
end

cdnaDb = [params(1).species '.cdnaDb.fas'];
% MatLab's use of blastlocal requires short entry names
simpleHeader = Header;
for n = 1:length(Header)
    pos = regexp(Header{n,1}, ':');
    simpleHeader{n,1} = Header{n,1}(1:pos(1)-1);
end
if exist(cdnaDb, 'file')
    delete([cdnaDb '*']);
end
fastawrite(cdnaDb, simpleHeader, Sequence);
blastformat('Inputdb', cdnaDb,...
    'FormatArgs', '-o T -p F');

