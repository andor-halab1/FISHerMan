function probeList=generateProbeList(adapterList,probeHeader,probeSequence,varargin)

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,'number',48);
end

if params.verbose
    disp('generating the list of probes to order');
    disp('  removing transcripts without enough oligos');
end

pos = regexp(probeHeader, 'ENS\w*T\d*', 'end');
trimHeader = {};
for n = 1:length(probeHeader)
    trimHeader{end+1} = probeHeader{n,1}(1:pos{n,1});
end
trimHeader = trimHeader';
uniqueHeader = unique(trimHeader, 'stable');

indexTotal = zeros(length(trimHeader),1);
for n = 1:length(uniqueHeader)
    index = ismember(trimHeader, uniqueHeader{n,1});
    if sum(index) < params.number && isempty(strfind(uniqueHeader{n,1},'ENSSPT')) ...
            && isempty(strfind(uniqueNames{n,1},'ENSMUST00000100497')) ...
            && isempty(strfind(uniqueNames{n,1},'ENSMUST00000118875')) % for Bin's special sequences
        indexTotal = indexTotal+index;
        if params.verbose
            disp(['  transcript ' uniqueHeader{n,1} ...
                ' has less than ' num2str(params.number) ' probes']);
        end
    end
end

indexTotal = logical(indexTotal);
probeHeader(indexTotal) = [];
probeSequence(indexTotal) = [];

disp('  randomizing and saving the list of probes');
index = randperm(length(probeHeader))';
probeHeader = probeHeader(index);
probeSequence = probeSequence(index);
probeList = [params.species '.probes.fas'];
if exist(probeList, 'file')
    delete(probeList);
end
fastawrite(probeList,probeHeader,probeSequence);

%% Check how many transcripts are left after this step of screening
[geneNumLeft,geneNumDelete]=checkTranscriptsLeft(adapterList,probeHeader);
if params.verbose
    disp([num2str(geneNumDelete) ' out of ' num2str(geneNumLeft+geneNumDelete)...
        ' FISH escaped FISHerMan''s net']);
end
