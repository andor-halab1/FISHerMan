function filePathList = blastFileSplit(Header, Sequence, params)

% params = struct('species','Mouse','verbose',1,'seqNum',1000);

seqNum = params(1).seqNum;
filePathList = {};

% MatLab's use of blastlocal requires short entry names
for n = 1:length(Header)
    pos = regexp(Header{n,1}, ':');
    if ~isempty(pos)
        Header{n,1} = Header{n,1}(1:pos(1)-1);
    else
        Header{n,1} = Header{n,1};
    end
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

