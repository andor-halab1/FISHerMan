function filePathList = blastFileSplit(Header, Sequence, seqNum, params)

% switch length(varargin)
%     case 0
%         seqNum = 48;
%         params = struct('species','Mouse','verbose',1,'keys','ENS\w*T\d*');
%     case 1
%         seqNum = varargin{1};
%         params = struct('species','Mouse','verbose',1,'keys','ENS\w*T\d*');
%     otherwise
%         seqNum = varargin{1};
%         params = varargin{2};
% end

filePathList = {};

% MatLab's use of blastlocal requires short entry names
for n = 1:length(Header)
    pos = regexp(Header{n,1}, params(1).keys, 'end');
    Header{n,1} = Header{n,1}(1:pos);
    Header{n,1} = strcat(Header{n,1}, '=', num2str(n));
end

delete('temp_*.fas');

for n = 1:floor(length(Header)/seqNum)
    tempName = ['temp_' num2str(n) '.fas'];
    filePathList{end+1} = tempName;
    fastawrite(tempName, Header((n-1)*seqNum+1:n*seqNum), Sequence((n-1)*seqNum+1:n*seqNum));
end

if floor(length(Header)/seqNum) < ceil(length(Header)/seqNum)
    n = ceil(length(Header)/seqNum);
    tempName = ['temp_' num2str(n) '.fas'];
    filePathList{end+1} = tempName;
    fastawrite(tempName, Header((n-1)*seqNum+1:end), Sequence((n-1)*seqNum+1:end));
end

filePathList = filePathList';

