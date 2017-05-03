function oligos = oligosParse(varargin)

% oligos = 'C:\OligoArray\oligos.txt';

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,'length',30,'number',48);
end

oligos = [params.species '.tempoligos.txt'];
if ~exist(oligos, 'file')
    warning('missing important files from OligoArray');
end

if params.verbose
    disp('reading the result file from OligoArray');
end

fid = fopen(oligos,'r');
fmt = '%s %f %f %f %f %f %f %s %s %s %*[^\n]';
temp = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','TreatAsEmpty','NA');
fclose(fid);

%% Remove non-specific oligos
if params.verbose
    disp('removing non-specific oligos');
end

geneNames = temp{1,1};
nonspecificHits = temp{1,3}(:,1);
specificHits = temp{1,3}(:,2);
index = [];
for n = 1:length(geneNames)
    if params.verbose && mod(n, 1000) == 1
        disp(['  analyzing oligo entry no. ' num2str(n)]);
    end
    [pos1, pos2] = regexp(nonspecificHits{n,1}, 'ENS\w*G\d*', 'start', 'end');
    
    flag = 0;
    for m = 1:length(pos1)
        if ~strfind(geneNames{n,1}, nonspecificHits{n,1}(pos1(m):pos2(m)))
            flag = 1;
        end
    end
    if flag == 1
        index = [index n];
    end
end
index = unique(index','stable');
% geneNames2 = geneNames;
% nonspecificHits2 = nonspecificHits;
% specificHits2 = specificHits;
geneNames(index) = [];
nonspecificHits(index) = [];
specificHits(index) = [];

%% Convert oligos into their complementary sequences
% OligoArray generates same sequences as transcriptome sequences, not
% their complementary
for n = 1:length(specificHits)
    specificHits{n,1} = seqrcomplement(specificHits{n,1});
end

%% Remove transcripts without enough oligos
% if params.verbose
%     disp('removing transcripts without enough oligos');
% end

% pos = regexp(geneNames, 'ENS\w*T\d*', 'end');
% trimNames = {};
% for n = 1:length(geneNames)
%     trimNames{end+1} = geneNames{n,1}(1:pos{n,1});
% end
% trimNames = trimNames';
% uniqueNames = unique(trimNames, 'stable');
% 
% indexTotal = zeros(length(trimNames),1);
% for n = 1:length(uniqueNames)
%     index = ismember(trimNames, uniqueNames{n,1});
%     if sum(index) < params.number
%         indexTotal = indexTotal+index;
%         disp(['transcript ' uniqueNames{n,1} ' has less than ' num2str(params.number) ' probes']);
%     end
% end
% 
% indexTotal = logical(indexTotal);
% geneNames(indexTotal) = [];
% nonspecificHits(indexTotal) = [];
% specificHits(indexTotal) = [];

%% Blast oligos against abundant rna database and remove non-specific oligos
[geneNames,specificHits,nonspecificHits]=blastAbundantRNASimple(geneNames,specificHits,nonspecificHits);

%% Remove transcripts without enough oligos
if params.verbose
    disp(['there are ' num2str(length(unique(geneNames))) ' transcripts']);
    disp('removing transcripts without enough oligos');
end

pos = regexp(geneNames, 'ENS\w*T\d*', 'end');
trimNames = {};
for n = 1:length(geneNames)
    trimNames{end+1} = geneNames{n,1}(1:pos{n,1});
end
trimNames = trimNames';
uniqueNames = unique(trimNames, 'stable');

indexTotal = zeros(length(trimNames),1);
for n = 1:length(uniqueNames)
    index = ismember(trimNames, uniqueNames{n,1});
    if sum(index) < params.number
        indexTotal = indexTotal+index;
        if params.verbose
            disp(['  transcript ' uniqueNames{n,1} ...
                ' has less than ' num2str(params.number) ' probes']);
        end
    end
end

indexTotal = logical(indexTotal);
geneNames(indexTotal) = [];
nonspecificHits(indexTotal) = [];
specificHits(indexTotal) = [];

%% Rearrange the list of oligos by grouping same genes together
if params.verbose
    disp('rearranging the oligo fasta file');
end

[geneNames,specificHits,nonspecificHits]=rearrangeOligos(geneNames,specificHits,nonspecificHits);

%% Save the list of oligos into files
if params.verbose
    disp('saving the oligo fasta file');
end

delete(oligos);
oligos = [params.species '.oligos.fas'];
if exist(oligos, 'file')
    delete(oligos);
end
fastawrite(oligos,geneNames,specificHits);

if params.verbose
    disp('done parsing OligoArray results');
end

