function [geneNumLeft,geneNumDelete]=checkTranscriptsLeft(adapterList,probeHeader,params)

% if length(varargin) >= 1
%     params = varargin{1};
% else
%     params = struct('species','Mouse','verbose',1,'keys','ENS\w*T\d*');
% end

[adapterHeader, adapterSequence] = fastaread(adapterList);
adapterHeader = adapterHeader';
adapterSequence = adapterSequence';
geneNumTotal = length(adapterHeader);

simpleHeader = probeHeader;
for n = 1:length(probeHeader)
    pos = regexp(probeHeader{n,1}, params(1).keys, 'end');
    simpleHeader{n,1} = probeHeader{n,1}(1:pos);
end
uniqueHeader = unique(simpleHeader,'stable');
geneNumLeft = length(uniqueHeader);
geneNumDelete = geneNumTotal-geneNumLeft;

[adapterHeader,adapterSequence] = pickExpressedSeq(uniqueHeader,adapterHeader,adapterSequence,params);
if exist(adapterList, 'file')
    delete(adapterList);
end
fastawrite(adapterList, adapterHeader, adapterSequence);

