function [ncrna,trna,Header,Sequence]=ncrnaParse(ncrna,trna,varargin)

% ncrna = 'C:\OligoArray\Mouse38.ncrna.fa';

switch length(varargin)
    case 0
        seqData = [];
        params = struct('species','Mouse','verbose',1,...
            'keys',{'ENS\w*T\d*','ENS\w*G\d*','gene_symbol:\S*','gene_biotype:\S*'});
    case 1
        seqData = varargin{1};
        params = struct('species','Mouse','verbose',1,...
            'keys',{'ENS\w*T\d*','ENS\w*G\d*','gene_symbol:\S*','gene_biotype:\S*'});
    otherwise
        seqData = varargin{1};
        params = varargin{2};
end

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
    temp1 = temp(pos1{1,1}{n,1}:pos2{1,1}{n,1});
    temp2 = temp(pos1{2,1}{n,1}:pos2{2,1}{n,1});
    temp3 = temp(pos1{3,1}{n,1}+12:pos2{3,1}{n,1});
    temp4 = temp(pos1{4,1}{n,1}+13:pos2{4,1}{n,1});
    
    if isempty(temp1)
        disp('missing transcript ID');
    end
    if isempty(temp2)
        disp('missing gene ID');
    end
    if isempty(temp3)
        disp('missing gene name');
    end
    if isempty(temp4)
        disp('missing gene type');
    end
    
    if isempty(strfind(temp,'pseudogene'))
        Header{n,1} = strcat(temp1, '-', temp2, '-', temp3, '-', temp4);
    else
        Header{n,1} = strcat(temp1, '-', temp2, '-', temp3, '-', temp4, '=pseudogene');
    end
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
        Header{end+1,1}=['ENSTRNAT00000000000-ENSTRNAG00000000000-tRNA' num2str(n) '-tRNA'];
        Sequence{end+1,1}=trnaSequence{n,1};
    end
    temp = strfind(trna,'\');
    trna = trna(temp(end)+1:end);
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
    pos = regexp(Header{n,1}, params(1).keys, 'end');
    simpleHeader{n,1} = Header{n,1}(1:pos);
end
if exist(ncrnaDb, 'file')
    delete([ncrnaDb '*']);
end
fastawrite(ncrnaDb, simpleHeader, Sequence);
blastformat('Inputdb', ncrnaDb,...
    'FormatArgs', '-o T -p F');

