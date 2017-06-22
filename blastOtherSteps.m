% function [probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
%     =blastOtherSteps(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore,params)
function [probeHeader,probeSequence,probeSequence2Seg,probeSequenceCore]...
    =blastOtherSteps(adapterList,probeHeader,probeSequence,probeSequence2Seg,probeSequenceCore,params)

% params = struct('species','Mouse','verbose',1,'keys',...
%     'thres',22,'querySize',50,...
%     'blastArgs1','-S 2','blastArgs2','-S 3',...
%     'grr','CCGCAACATCCAGCATCGTG',...
%     'rRr','AGAGTGAGTAGTAGTGGAGT','rGr','GATGATGTAGTAGTAAGGGT',...
%     'rBr','TGTGATGGAAGTTAGAGGGT','rIRr','GGAGTAGTTGGTTGTTAGGA');

if params(1).verbose
    disp('removing probes that non-specifically bind to 2nd PCR primers and other probes');
end

[adapterHeader, adapterSequence] = fastaread(adapterList);
adapterHeader = adapterHeader';
adapterSequence = adapterSequence';

simpleHeader = probeHeader;
for n = 1:length(probeHeader)
    pos = regexp(probeHeader{n,1}, ':');
    simpleHeader{n,1} = probeHeader{n,1}(1:pos(1)-1);
end

seqDelete = [];
for n = 1:length(adapterHeader)
    if params(1).verbose% && mod(n, 1000) == 1
        disp(['  working on trancript ' adapterHeader{n,1}]);
    end
    temp...
        =blastOneTranscript(adapterHeader{n,1},adapterSequence{n,1},simpleHeader,probeSequenceCore,params);
    if ~isempty(temp)
        seqDelete = [seqDelete temp];
    end
end

probeHeader(seqDelete)= [];
probeSequence(seqDelete)= [];
% probeSequence3Seg(seqDelete)= [];
probeSequence2Seg(seqDelete)= [];
probeSequenceCore(seqDelete)= [];

%% Check how many transcripts are left after this step of screening
[geneNumLeft,geneNumDelete]=checkTranscriptsLeft(adapterList,probeHeader);
if params(1).verbose
    disp([num2str(geneNumDelete) ' out of ' num2str(geneNumLeft+geneNumDelete)...
        ' FISH escaped FISHerMan''s net']);
end

