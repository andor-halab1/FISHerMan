function [geneNumLeft,geneNumDelete]=checkTranscriptsLeft(adapterList,probeHeader)

% params = struct('species','Mouse','verbose',1);

[adapterHeader, adapterSequence] = fastaread(adapterList);
adapterHeader = adapterHeader';
adapterSequence = adapterSequence';
geneNumTotal = length(adapterHeader);

simpleHeader = probeHeader;
for n = 1:length(probeHeader)
    pos = regexp(probeHeader{n,1}, ':');
    simpleHeader{n,1} = probeHeader{n,1}(1:pos(1)-1);
end
uniqueHeader = unique(simpleHeader,'stable');
geneNumLeft = length(uniqueHeader);
geneNumDelete = geneNumTotal-geneNumLeft;

[adapterHeader,adapterSequence] = pickExpressedSeq(uniqueHeader,adapterHeader,adapterSequence);
if exist(adapterList, 'file')
    delete(adapterList);
end
fastawrite(adapterList, adapterHeader, adapterSequence);

