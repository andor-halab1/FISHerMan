function [probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =blast1stPCR(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore,varargin)

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,'keys','ENS\w*T\d*',...
        'thres',22,'querySize',20,'seqNum',1000,...
        'blastArgs','','parallel', 0,...
        'gf','GGAATCGTTGCGGGTGTCCT','grr','CCGCAACATCCAGCATCGTG',...
        'T7r','CCCTATAGTGAGTCGTATTA',...
        'rRr','AGAGTGAGTAGTAGTGGAGT','rGr','GATGATGTAGTAGTAAGGGT',...
        'rBr','TGTGATGGAAGTTAGAGGGT','rIRr','GGAGTAGTTGGTTGTTAGGA');
end

if params.verbose
    disp('removing probes that non-specifically bind to primers in the 1st PCR step');
end

%% Generate probe database files for Blast
if params.verbose
    disp('  generating probe database files for Blast');
end

probesDb = [params.species '.probesDb.fas'];
% MatLab's use of blastlocal requires short entry names
simpleHeader = probeHeader;
for n = 1:length(probeHeader)
    pos = regexp(probeHeader{n,1}, params.keys, 'end');
    simpleHeader{n,1} = strcat(probeHeader{n,1}(1:pos),'=',num2str(n));
end
if exist(probesDb, 'file')
    delete([probesDb '*']);
end
fastawrite(probesDb, simpleHeader, probeSequenceCore);
blastformat('Inputdb', probesDb,...
    'FormatArgs', '-o T -p F');

params.DbSize = getDbSize(probesDb);

%% Split one giant fasta file into smaller ones, so that parallel computing is possible
if params.verbose
    disp('  spliting fasta files for parallel computing');
end

[adapterHeader, adapterSequence] = fastaread(adapterList);
adapterHeader = adapterHeader';
adapterSequence = adapterSequence';
adapterHeader{end+1,1} = 'ENSPRIMERT00';
adapterSequence{end+1,1} = params.gf;
adapterHeader{end+1,1} = 'ENSPRIMERT01';
adapterSequence{end+1,1} = params.grr;

filePathList = blastFileSplit(adapterHeader, adapterSequence, params.seqNum);
fileNum = length(filePathList);

%% Blast primers against probes
eValue = bitScore2eValue(params.thres, params.querySize, params.DbSize);

DbPath = probesDb;
blastArgs = [params.blastArgs ' -e ' num2str(eValue) ' -b ' num2str(length(probeHeader))];

blastData = {};
if params.parallel
    poolobj = parpool;
    verbose = params.verbose;
    parfor k = 1:fileNum
        if verbose
            disp(['  blasting temporary file no. ' num2str(k)]);
            startTime(k) = tic;
        end
        blastData{k,1} = blastOp(filePathList{k}, DbPath, blastArgs);
        if verbose
            totalTime(k) = toc(startTime(k));
            disp(['  elapsed time is ' num2str(totalTime(k)) ' seconds']);
        end
    end
    delete(poolobj);
else
    for k = 1:fileNum
        if params.verbose
            disp(['  blasting temporary file no. ' num2str(k)]);
            startTime = tic;
        end
        blastData{k,1} = blastOp(filePathList{k}, DbPath, blastArgs);
        if params.verbose
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
        seqDelete = [seqDelete probeNum];
    end
end
seqDelete = seqDelete';

probeHeader(seqDelete)= [];
probeSequence(seqDelete)= [];
probeSequence3Seg(seqDelete)= [];
probeSequenceCore(seqDelete)= [];

%% Check how many transcripts are left after this step of screening
[geneNumLeft,geneNumDelete]=checkTranscriptsLeft(adapterList,probeHeader);
if params.verbose
    disp([num2str(geneNumDelete) ' out of ' num2str(geneNumLeft+geneNumDelete)...
        ' FISH escaped FISHerMan''s net']);
end

fileCleanUp(filePathList);

