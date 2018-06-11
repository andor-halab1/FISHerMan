function oligos = oligosParse(params)

% params = struct('species','Mouse','verbose',1,...
%     'number',48,'seqNum',1000,'thres',30,'querySize',30,...
%     'blastArgs','-S 2','parallel', 0,...
%     'specialTranscripts','C:\FISHerMan\Db\Mouse.STList.fas');

if params(1).verbose
    disp('Reading the result file from OligoArray');
end

oligos = [params(1).species '.tempoligos.txt'];
if exist(oligos, 'file')
    fid = fopen(oligos,'r');
    fmt = '%s %f %f %f %f %f %f %s %s %s %*[^\n]';
    temp = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','TreatAsEmpty','NA');
    fclose(fid);
else
    warning('Missing important files from OligoArray');
    return;
end

% Remove non-specific oligos
if params(1).verbose
    disp('Removing non-specific oligos');
end

geneNames = temp{1,1};
nonspecificHits = temp{1,3}(:,1);
specificHits = temp{1,3}(:,2);
index = [];
for n = 1:length(geneNames)
    if params(1).verbose && mod(n, 1000) == 1
        disp(['  analyzing oligo entry no. ' num2str(n)]);
    end
    
    % all Header here are two-segment Header
    pos = regexp(geneNames{n,1}, ':');
    if ~isempty(pos)
        geneName=geneNames{n,1}(pos(1)+1:end);
    end
    
    if length(regexp(nonspecificHits{n,1}, geneName)) < ...
            length(regexp(nonspecificHits{n,1}, ':'))
        index = [index;n];
    end
end

index = unique(index); 
geneNames(index) = [];
nonspecificHits(index) = [];
specificHits(index) = [];

%% Convert oligos into their complementary sequences
% OligoArray generates same sequences as transcriptome sequences, not
% their complementary
for n = 1:length(specificHits)
    specificHits{n,1} = seqrcomplement(specificHits{n,1});
end

%% Blast oligos against abundant rna database and remove non-specific oligos
[geneNames,specificHits,nonspecificHits,~] = blastAbundantRNA([], geneNames, specificHits, nonspecificHits, [], params);


%% Remove transcripts without enough oligos
if params(1).verbose
    disp(['there are ' num2str(length(unique(geneNames))) ' transcripts']);
    disp('removing transcripts without enough oligos');
end

pos = regexp(geneNames, ':');
trimNames = {};
for n = 1:length(geneNames)
    trimNames{end+1,1} = geneNames{n,1}(1:pos{n,1}(1)-1);
end
uniqueNames = unique(trimNames, 'stable');

indexTotal = zeros(length(trimNames),1);
for n = 1:length(uniqueNames)
    index = ismember(trimNames, uniqueNames{n,1});
    if sum(index) < params(1).number &&...
            checkSpecialTranscripts(uniqueNames{n,1},params) % check for Bin's special sequences
        indexTotal = indexTotal + index;
        if params(1).verbose
            disp(['  transcript ' uniqueNames{n,1} ' has less than ' num2str(params(1).number) ' probes']);
        end
    end
end

indexTotal = logical(indexTotal);
geneNames(indexTotal) = [];
nonspecificHits(indexTotal) = [];
specificHits(indexTotal) = [];

% Rearrange the list of oligos by grouping same genes together
if params(1).verbose
    disp('Rearranging the oligo fasta file');
end

[geneNames,specificHits,nonspecificHits] = rearrangeOligos(geneNames,specificHits,nonspecificHits,params);

% Save the list of oligos into files
if params(1).verbose
    disp('Saving the oligo fasta file');
end

oligos = [params(1).species '.oligos.txt'];
if exist(oligos, 'file')
    delete(oligos);
end
fastawrite(oligos,geneNames,specificHits);

if params(1).verbose
    disp('Done parsing OligoArray results');
end
