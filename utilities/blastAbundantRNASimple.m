% blast probes against abundant rna database. Only blast against the 
% complementary sequences of abundant rna. Rule out the possibility of
% homology of 15 nt or more.

function [Header,Sequence,nonSequence]...
    =blastAbundantRNASimple(Header,Sequence,nonSequence,params)

% params = struct('species','Mouse','verbose',1,...
%     'number',48,'seqNum',1000,'thres',30,'querySize',30,...
%     'DbSize',2*10^5,'blastArgs','-S 2','parallel', 0,...
%     'specialTranscripts','C:\FISHerMan\Db\Mouse.STList.fas');

if isempty(nonSequence)
    nonSequence = Sequence;
end

if params(1).verbose
    disp('removing non-specific oligos against the abundant rna database');
end

%% Split one giant fasta file into smaller ones, so that parallel computing is possible
if params(1).verbose
    disp('  spliting fasta files for parallel computing');
end

filePathList = blastFileSplit(Header, Sequence, params);
fileNum = length(filePathList);

%% Blast mouse oligos against abundant rna
DbPath = [params(1).species '.abundantrnaDb.fas'];
params(1).DbSize = getDbSize(DbPath);

eValue = bitScore2eValue(params(1).thres, params(1).querySize, params(1).DbSize);
blastArgs = [params(1).blastArgs ' -e ' num2str(eValue)];

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

%% Find out the oligos that non-specifically hit abundant rna database
seqDelete = [];
for n = 1:length(data)
    flag = 0;
    for m = 1:length(data(n).Hits)
        if isempty(strfind(data(n).Query,data(n).Hits(m).Name))
            flag = 1;
        end
    end
    if flag == 1
        seqDelete = [seqDelete n];
    end
end
seqDelete = seqDelete';

Header(seqDelete)= [];
Sequence(seqDelete)= [];
nonSequence(seqDelete)= [];

fileCleanUp(filePathList);

