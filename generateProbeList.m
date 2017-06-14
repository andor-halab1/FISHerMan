function probeList=generateProbeList(adapterList,probeHeader,probeSequence,params)

% params = struct('species','Mouse','verbose',1,'number',48,...
%     'specialTranscripts','C:\FISHerMan\Db\Mouse.STList.fas');

%% Remove transcripts without enough probes
if params(1).verbose
    disp('generating the list of probes to order');
    disp('  removing transcripts without enough probes');
end

pos = regexp(probeHeader, ':');
trimHeader = {};
for n = 1:length(probeHeader)
    trimHeader{end+1} = probeHeader{n,1}(1:pos{n,1}(1)-1);
end
trimHeader = trimHeader';
uniqueHeader = unique(trimHeader, 'stable');

indexTotal = zeros(length(trimHeader),1);
for n = 1:length(uniqueHeader)
    index = ismember(trimHeader, uniqueHeader{n,1});
    if sum(index) < params(1).number &&...
            checkSpecialTranscripts(uniqueHeader{n,1},params) % for Bin's special sequences
        indexTotal = indexTotal+index;
        if params(1).verbose
            disp(['  transcript ' uniqueHeader{n,1} ...
                ' has less than ' num2str(params(1).number) ' probes']);
        end
    end
end

indexTotal = logical(indexTotal);
probeHeader(indexTotal) = [];
probeSequence(indexTotal) = [];
trimHeader(indexTotal) = [];

%% Trim the number of probes of each transcript to be 48
uniqueHeader = unique(trimHeader, 'stable');
indexTotal = [];
for n = 1:length(uniqueHeader)
    index = ismember(trimHeader, uniqueHeader{n,1});
    if sum(index) > params(1).number
        index=find(index>0);
        index=index(params(1).number+1:end);
        indexTotal=[indexTotal;index];
    end
end

probeHeader(indexTotal) = [];
probeSequence(indexTotal) = [];

% disp('  randomizing and saving the list of probes');
% indexTotal = randperm(length(probeHeader))';
% probeHeader = probeHeader(indexTotal);
% probeSequence = probeSequence(indexTotal);
probeList = [params(1).species '.probes.nr.txt'];
if exist(probeList, 'file')
    delete(probeList);
end
fastawrite(probeList,probeHeader,probeSequence);

%% Check how many transcripts are left after this step of screening
[geneNumLeft,geneNumDelete]=checkTranscriptsLeft(adapterList,probeHeader);
if params(1).verbose
    disp([num2str(geneNumDelete) ' out of ' num2str(geneNumLeft+geneNumDelete)...
        ' FISH escaped FISHerMan''s net']);
end
