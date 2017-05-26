function [expressedHeader, expressedSequence] = pickExpressedSeq(seqData, Header, Sequence, params)

% if length(varargin) >= 1
%     params = varargin{1};
% else
%     params = struct('species','Mouse','verbose',1,'keys','ENS\w*T\d*');
% end

transcriptID = seqData(:,1);

pos = regexp(Header, params(1).keys, 'end');
temp = Header;
for n = 1:length(Header)
    temp{n,1} = Header{n,1}(1:pos{n,1});
end

[~,index,~] = intersect(temp, transcriptID, 'stable');
expressedHeader = Header(index);
expressedSequence = Sequence(index);
