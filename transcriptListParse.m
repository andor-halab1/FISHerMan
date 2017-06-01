function [longHeader,longSequence]=transcriptListParse(transcriptList,cdnaHeader,cdnaSequence,...
    ncrnaHeader,ncrnaSequence,params)

% params = struct('species','Mouse','verbose',1,...
%     'dir1','C:\FISHerMan\Db\Mouse.transcriptList.fas',...
%     'length',40,'number',24);

%% Pick transcript sequences that are longer than a certain threshold
if params(1).verbose
    disp('reading the name file of target transcripts');
    disp('  thresholding based on sequence length');
end

% use provided transcriptList as targets for designing FISH probes
[targetHeader,~]=fastaread(transcriptList);
targetHeader = targetHeader';

totalHeader = vertcat(cdnaHeader,ncrnaHeader);
totalSequence = vertcat(cdnaSequence,ncrnaSequence);
[totalHeader, totalSequence] = pickExpressedSeq(targetHeader, totalHeader, totalSequence);

longHeader = {};
longSequence = {};
for n = 1:length(totalSequence)
    if length(totalSequence{n,1}) >= (params(1).length*params(1).number)
        if length(regexp(totalHeader{n,1},':')) >= 2
            pos=regexp(totalHeader{n,1},':');
            totalHeader{n,1}=totalHeader{n,1}(1:pos(2)-1);
        end
        longHeader{end+1,1} = totalHeader{n,1};
        longSequence{end+1,1} = totalSequence{n,1};
    end
end

if params(1).verbose
    disp(['  ' num2str(length(longHeader)) ' transcripts have length longer than '...
        num2str(params(1).length*params(1).number) ' nt']);
end

%% Segment each transcript sequence into 1kb pieces, so later analysis runs more efficiently
if params(1).verbose
    disp('  segmenting long sequences into 1kb segments');
end

segHeader = {};
segSequence = {};
for n = 1:length(longHeader)
    segN = floor(length(longSequence{n,1})/1000);
    probe = longSequence{n,1};
    for m = 1:segN
        segHeader{end+1,1} = longHeader{n,1};
        segSequence{end+1,1} = probe(1000*(m-1)+1:1000*m);
    end
    if length(probe) >= (1000*segN+30)
        segHeader{end+1,1} = longHeader{n,1};
        segSequence{end+1,1} = probe(1000*segN+1:end);
    end
end

if params(1).verbose
    disp('saving the target transcript fasta file');
end

transcriptList = [params(1).species '.target1kb.fas'];
if exist(transcriptList, 'file')
    delete(transcriptList);
end
fastawrite(transcriptList, segHeader, segSequence);

