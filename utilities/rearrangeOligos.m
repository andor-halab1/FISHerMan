% Rearrange the list of oligos by grouping same genes together.

function [rearrangedHeader,rearrangedSequence,rearrangednonSequence]=...
    rearrangeOligos(Header,Sequence,nonSequence,params)

% if length(varargin) >= 1
%     params = varargin{1};
% else
%     params = struct('species','Mouse','verbose',1);
% end

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

