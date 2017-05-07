function [transcriptList,longHeader,longSequence]=transcriptListParse(transcriptList,cdnaHeader,cdnaSequence,...
    ncrnaHeader,ncrnaSequence,varargin)
%% Pick transcript sequences that are longer than a certain threshold
% transcriptList = 'C:\FISHerMan\transcriptList.fas';

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,'length',40,'number',24);
end

% if params.verbose
%     disp('generating files that contain the list of target transcripts');
% end

if ~isempty(transcriptList) % use provided seqList as targets for designing FISH probes
    if params.verbose
        disp('reading the name file of target transcripts');
        disp('  thresholding based on sequence length');
    end

    [targetHeader,~]=fastaread(transcriptList);
    targetHeader = targetHeader';
    
    longHeader = {};
    longSequence = {};
    for n = 1:length(cdnaSequence)
        temp = cdnaHeader{n,1};
        pos = regexp(temp, 'ENS\w*T\d*', 'end');
        temp = temp(1:pos);
        if length(cdnaSequence{n,1}) >= (params.length*params.number) &&...
                sum(ismember(targetHeader, temp))
            longHeader{end+1} = cdnaHeader{n,1};
            longSequence{end+1} = cdnaSequence{n,1};
        end
    end
    for n = 1:length(ncrnaSequence)
        temp = ncrnaHeader{n,1};
        pos = regexp(temp, 'ENS\w*T\d*', 'end');
        temp = temp(1:pos);
        if length(ncrnaSequence{n,1}) >= (params.length*params.number) &&...
                sum(ismember(targetHeader, temp))
            longHeader{end+1} = ncrnaHeader{n,1};
            longSequence{end+1} = ncrnaSequence{n,1};
        end
    end
else % use all possible transcripts in the cdna file as targets for designing FISH probes
    if params.verbose
        disp('using all possible transcripts in the cdna fasta file');
        disp('  thresholding based on sequence length');
    end
    
    pos =  regexp(cdnaHeader, 'pseudogene');
    
    longHeader = {};
    longSequence = {};
    for n = 1:length(cdnaHeader)
        if length(cdnaSequence{n,1}) >= (params.length*params.number) && isempty(pos{n,1})
            longHeader{end+1} = cdnaHeader{n,1};
            longSequence{end+1} = cdnaSequence{n,1};
        end
    end
end
longHeader = longHeader';
longSequence = longSequence';

if params.verbose
    disp(['  ' num2str(length(longHeader)) ' transcripts have length longer than '...
        num2str(params.length*params.number) ' nt']);
end

%% Segment each transcript sequence into 1kb pieces, so later analysis runs more efficiently
if params.verbose
    disp('  segmenting long sequences into 1kb segments');
end

segHeader = {};
segSequence = {};
for n = 1:length(longHeader)
    segN = floor(length(longSequence{n,1})/1000);
    probe = longSequence{n,1};
    for m = 1:segN
        segHeader{end+1} = longHeader{n,1};
        segSequence{end+1} = probe(1000*(m-1)+1:1000*m);
    end
    if length(probe) >= (1000*segN+30)
        segHeader{end+1} = longHeader{n,1};
        segSequence{end+1} = probe(1000*segN+1:end);
    end
end
segHeader = segHeader';
segSequence = segSequence';

if params.verbose
    disp('saving the target transcript fasta file');
end

transcriptList = [params(1).species '.target1kb.fas'];
if exist(transcriptList, 'file')
    delete(transcriptList);
end
fastawrite(transcriptList, segHeader, segSequence);

