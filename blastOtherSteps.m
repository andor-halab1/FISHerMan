function [probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =blastOtherSteps(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore,varargin)

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,'keys','ENS\w*T\d*',...
        'thres',22,'querySize',50,'blastArgs','',...
        'grr','CCGCAACATCCAGCATCGTG','T7r','CCCTATAGTGAGTCGTATTA',...
        'rRr','AGAGTGAGTAGTAGTGGAGT','rGr','GATGATGTAGTAGTAAGGGT',...
        'rBr','TGTGATGGAAGTTAGAGGGT','rIRr','GGAGTAGTTGGTTGTTAGGA');
end

if params.verbose
    disp('removing probes that non-specifically bind to 2nd PCR primers and other probes');
end

[adapterHeader, adapterSequence] = fastaread(adapterList);
adapterHeader = adapterHeader';
adapterSequence = adapterSequence';

simpleHeader = probeHeader;
for n = 1:length(probeHeader)
    pos = regexp(probeHeader{n,1}, params.keys, 'end');
    simpleHeader{n,1} = probeHeader{n,1}(1:pos);
end

seqDelete = [];
for n = 1:length(adapterHeader)
    if params.verbose% && mod(n, 1000) == 1
        disp(['  working on trancript ' adapterHeader{n,1}]);
    end
    temp...
        =blastOneTranscript(adapterHeader{n,1},adapterSequence{n,1},simpleHeader,probeSequenceCore);
    seqDelete = [seqDelete temp];
end

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

