function seqDelete...
    =blastOneTranscript(OTAdapterHeader,OTAdapterSequence,probeHeader,probeSequenceCore,varargin)

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,...
        'thres',22,'querySize',50,'blastArgs','',...
        'grr','CCGCAACATCCAGCATCGTG','T7r','CCCTATAGTGAGTCGTATTA',...
        'rRr','AGAGTGAGTAGTAGTGGAGT','rGr','GATGATGTAGTAGTAAGGGT',...
        'rBr','TGTGATGGAAGTTAGAGGGT','rIRr','GGAGTAGTTGGTTGTTAGGA');
end

index = ismember(probeHeader, OTAdapterHeader);
temp = 1:length(probeHeader);
index = temp(index);
OTHeader = probeHeader(index);
OTSequenceCore = probeSequenceCore(index);

%% Generate database files for Blast
OTDb1 = [params.species '.' OTAdapterHeader '.Db1.fas'];
OTDb2 = [params.species '.' OTAdapterHeader '.Db2.fas'];
% MatLab's use of blastlocal requires short entry names

for n = 1:length(OTHeader)
    OTHeader{n,1} = strcat(OTHeader{n,1},'=',num2str(n));
end

Header{1,1} = 'ENSPRIMERT10';
Sequence{1,1} = strcat(params.rRr,OTAdapterSequence);
Header{2,1} = 'ENSPRIMERT11';
Sequence{2,1} = strcat(params.rGr,OTAdapterSequence);
Header{3,1} = 'ENSPRIMERT12';
Sequence{3,1} = strcat(params.rBr,OTAdapterSequence);
Header{4,1} = 'ENSPRIMERT13';
Sequence{4,1} = strcat(params.rIRr,OTAdapterSequence);
Header{5,1} = 'ENSPRIMERT14';
Sequence{5,1} = strcat(params.grr,params.T7r);

if exist(OTDb1, 'file')
    delete([OTDb1 '*']);
end
fastawrite(OTDb1, OTHeader, OTSequenceCore);
blastformat('Inputdb', OTDb1,...
    'FormatArgs', '-o T -p F');

if exist(OTDb2, 'file')
    delete([OTDb2 '*']);
end
fastawrite(OTDb2, Header, Sequence);
blastformat('Inputdb', OTDb2,...
    'FormatArgs', '-o T -p F');

params.DbSize1 = getDbSize(OTDb1);
params.DbSize2 = getDbSize(OTDb2);

%% Generate a temporary probe list file for blast
filePathList = blastFileSplit(OTHeader, OTSequenceCore, length(OTHeader));

%% Blast probes against 2nd PCR primers and other probes
eValue1 = bitScore2eValue(params.thres, params.querySize, params.DbSize1);
eValue2 = bitScore2eValue(params.thres, params.querySize, params.DbSize2);

DbPath1 = OTDb1;
blastArgs1 = [params.blastArgs ' -e ' num2str(eValue1) ' -b ' num2str(length(OTHeader)) ' -S 2'];

DbPath2 = OTDb2;
blastArgs2 = [params.blastArgs ' -e ' num2str(eValue2) ' -b ' num2str(length(OTHeader))];

if params.verbose
    disp('  blasting the temporary probe list file');
    startTime = tic;
end
blastData{1,1} = blastOp(filePathList{1}, DbPath1, blastArgs1);
blastData{2,1} = blastOp(filePathList{1}, DbPath2, blastArgs2);
if params.verbose
    totalTime = toc(startTime);
    disp(['  elapsed time is ' num2str(totalTime) ' seconds']);
end

%% Find out the probes that non-specifically bind 2nd PCR primers and other probes
seqDelete = findSeqDelete(blastData);

seqDelete = index(seqDelete);

fileCleanUp(filePathList);
delete([OTDb1 '*']);
delete([OTDb2 '*']);

