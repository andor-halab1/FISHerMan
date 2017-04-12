% blast probes against abundant rna database. Only blast against the 
% complementary sequences of abundant rna. Rule out the possibility of
% homology of 15 nt or more.

function [Header,Sequence,nonSequence,nonSequence2]...
    =blastAbundantRNA(adapterList,Header,Sequence,nonSequence,nonSequence2,varargin)

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,...
        'thres',30,'querySize',73,'DbSize',2*10^5,'seqNum',1000,...
        'blastArgs','-S 2','parallel', 0);
end

if isempty(nonSequence)
    nonSequence = Sequence;
end

if isempty(nonSequence2)
    nonSequence2 = Sequence;
end

if params.verbose
    disp('removing non-specific probes against the abundant rna database');
end

%% Split one giant fasta file into smaller ones, so that parallel computing is possible
if params.verbose
    disp('  spliting fasta files for parallel computing');
end

filePathList = blastFileSplit(Header, Sequence, params.seqNum);
fileNum = length(filePathList);

%% Blast mouse oligos against abundant rna
eValue = bitScore2eValue(params.thres, params.querySize, params.DbSize);

DbPath = [params.species '.abundantrnaDb.fas'];
blastArgs = [params.blastArgs ' -e ' num2str(eValue)];

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

%% Find out the oligos that non-specifically hit abundant rna database
seqDelete = [];
for n = 1:length(data)
    if ~isempty(data(n).Hits)
        seqDelete = [seqDelete n];
    end
end
seqDelete = seqDelete';

Header(seqDelete)= [];
Sequence(seqDelete)= [];
nonSequence(seqDelete)= [];
nonSequence2(seqDelete)= [];

%% Check how many transcripts are left after this step of screening
[geneNumLeft,geneNumDelete]=checkTranscriptsLeft(adapterList,Header);
if params.verbose
    disp([num2str(geneNumDelete) ' out of ' num2str(geneNumLeft+geneNumDelete)...
        ' FISH escaped FISHerMan''s net']);
end

fileCleanUp(filePathList);

