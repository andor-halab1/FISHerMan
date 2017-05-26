function [Header,Sequence]=cdnaParse(cdna,seqData,params)

% cdna = 'C:\OligoArray\Mouse38.cdna.fa';

% switch length(varargin)
%     case 0
%         seqData = [];
%         params = struct('species','Mouse','verbose',1,...
%             'keys',{'ENS\w*T\d*','ENS\w*G\d*'});
%     case 1
%         seqData = varargin{1};
%         params = struct('species','Mouse','verbose',1,...
%             'keys',{'ENS\w*T\d*','ENS\w*G\d*'});
%     otherwise
%         seqData = varargin{1};
%         params = varargin{2};
% end

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
    temp1 = temp(pos1{1,1}{n,1}:pos2{1,1}{n,1});
    temp2 = temp(pos1{2,1}{n,1}:pos2{2,1}{n,1});
    
    if isempty(temp1)
        disp('missing transcript ID');
    end
    if isempty(temp2)
        disp('missing gene ID');
    end
    Header{n,1} = strcat(temp1, ':', temp2);
end

if ~isempty(seqData)
    if params(1).verbose
        disp('  picking expressed sequences according to RNA-seq data');
    end
    [Header, Sequence] = pickExpressedSeq(seqData, Header, Sequence, params);
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
    pos = regexp(Header{n,1}, params(1).keys, 'end');
    simpleHeader{n,1} = Header{n,1}(1:pos);
end
if exist(cdnaDb, 'file')
    delete([cdnaDb '*']);
end
fastawrite(cdnaDb, simpleHeader, Sequence);
blastformat('Inputdb', cdnaDb,...
    'FormatArgs', '-o T -p F');

