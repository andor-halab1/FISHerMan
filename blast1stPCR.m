function [probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =blast1stPCR(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore,params)

% params = struct('species','Mouse','verbose',1,...
%     'thres',22,'querySize',20,'seqNum',1000,...
%     'blastArgs','-S 3','parallel', 0,...
%     'gf','GGAATCGTTGCGGGTGTCCT','grr','CCGCAACATCCAGCATCGTG');

if params(1).verbose
    disp('removing probes that non-specifically bind to primers in the 1st PCR step');
end

%% Generate probe database files for Blast
if params(1).verbose
    disp('  generating probe database files for Blast');
end

probesDb = [params(1).species '.probesDb.fas'];
% MatLab's use of blastlocal requires short entry names
simpleHeader = probeHeader;
for n = 1:length(probeHeader)
    pos = regexp(probeHeader{n,1}, ':');
    simpleHeader{n,1} = strcat(probeHeader{n,1}(1:pos(1)-1),'=',num2str(n));
end
if exist(probesDb, 'file')
    delete([probesDb '*']);
end
fastawrite(probesDb, simpleHeader, probeSequenceCore);
blastformat('Inputdb', probesDb,...
    'FormatArgs', '-o T -p F');

%% Split one giant fasta file into smaller ones, so that parallel computing is possible
if params(1).verbose
    disp('  spliting fasta files for parallel computing');
end

[adapterHeader, adapterSequence] = fastaread(adapterList);
adapterHeader = adapterHeader';
adapterSequence = adapterSequence';
adapterHeader{end+1,1} = 'ENSPRIMERT00';
adapterSequence{end+1,1} = params(1).gf;
adapterHeader{end+1,1} = 'ENSPRIMERT01';
adapterSequence{end+1,1} = params(1).grr;

filePathList = blastFileSplit(adapterHeader, adapterSequence, params);
fileNum = length(filePathList);

%% Blast primers against probes
DbPath = probesDb;
params(1).DbSize = getDbSize(DbPath);

eValue = bitScore2eValue(params(1).thres, params(1).querySize, params(1).DbSize);
blastArgs = [params(1).blastArgs ' -e ' num2str(eValue) ' -b ' num2str(length(probeHeader))];

blastData = {};
if params(1).parallel
    poolobj = parpool;
    verbose = params(1).verbose;
    parfor k = 1:fileNum
        if verbose
            disp(['  blasting temporary file no. ' num2str(k)]);
        end
        blastData{k,1} = blastOp(filePathList{k}, DbPath, blastArgs);
    end
    delete(poolobj);
else
    for k = 1:fileNum
        if params(1).verbose
            disp(['  blasting temporary file no. ' num2str(k)]);
            startTime = tic;
        end
        blastData{k,1} = blastOp(filePathList{k}, DbPath, blastArgs);
        if params(1).verbose
            totalTime = toc(startTime);
            disp(['  elapsed time is ' num2str(totalTime) ' seconds']);
        end
    end
end

data = [];
for k = 1:length(blastData)
    data = [data blastData{k,1}];
end

%% Find out the probes that non-specifically bind primers in the 1st PCR step
seqDelete = [];
for n = 1:length(data)
    for m = 1:length(data(n).Hits)
        probeName = data(n).Hits(m).Name;
        [probeNumPos1, probeNumPos2] = regexp(probeName, '=\d*', 'start', 'end');
        probeNum = probeName(probeNumPos1+1:probeNumPos2);
        probeNum = uint32(str2double(probeNum));
        seqDelete = [seqDelete;probeNum];
    end
end
seqDelete = unique(seqDelete);

probeHeader(seqDelete)= [];
probeSequence(seqDelete)= [];
probeSequence3Seg(seqDelete)= [];
probeSequenceCore(seqDelete)= [];

%% Check how many transcripts are left after this step of screening
[geneNumLeft,geneNumDelete]=checkTranscriptsLeft(adapterList,probeHeader);
if params(1).verbose
    disp([num2str(geneNumDelete) ' out of ' num2str(geneNumLeft+geneNumDelete)...
        ' FISH escaped FISHerMan''s net']);
end

fileCleanUp(filePathList);

