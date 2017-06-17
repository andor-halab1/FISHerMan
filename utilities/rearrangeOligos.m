% Rearrange the list of oligos by grouping same genes together.

function [rearrangedHeader,rearrangedSequence,rearrangednonSequence]=...
    rearrangeOligos(Header,Sequence,nonSequence,params)

% params = struct('species','Mouse','verbose',1,...
%     'number',48,'seqNum',1000,'thres',30,'querySize',30,...
%     'DbSize',2*10^5,'blastArgs','-S 2','parallel', 0,...
%     'specialTranscripts','C:\FISHerMan\Db\Mouse.STList.fas');

uniqueHeader = unique(Header, 'stable');

rearrangedHeader = {};
rearrangedSequence = {};
rearrangednonSequence = {};
for n = 1:length(uniqueHeader)
    if params(1).verbose && mod(n, 1000) == 1
        disp(['  rearranging transcript no. ' num2str(n)]);
    end
    index = ismember(Header, uniqueHeader{n,1});
    rearrangedHeader(end+1:end+sum(index),1) = Header(index);
    rearrangedSequence(end+1:end+sum(index),1) = Sequence(index);
    rearrangednonSequence(end+1:end+sum(index),1) = nonSequence(index);
end

