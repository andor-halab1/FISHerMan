function filePathList = blastFileSplit(Header, Sequence, params)

% params = struct('species','Mouse','verbose',1,...
%     'number',48,'seqNum',1000,'thres',30,'querySize',30,...
%     'DbSize',2*10^5,'blastArgs','-S 2','parallel', 0,...
%     'specialTranscripts','C:\FISHerMan\Db\Mouse.STList.fas');

seqNum = params(1).seqNum;
filePathList = {};

% MatLab's use of blastlocal requires short entry names
for n = 1:length(Header)
    pos1 = regexp(Header{n,1}, ':');
    pos2 = regexp(Header{n,1}, '=');
    if ~isempty(pos1)
        Header{n,1} = Header{n,1}(1:pos1(1)-1);
    elseif ~isempty(pos2)
        Header{n,1} = Header{n,1}(1:pos2(1)-1);
    end
    Header{n,1} = strcat(Header{n,1}, '=', num2str(n));
end

delete('temp_*.fas');

for n = 1:floor(length(Header)/seqNum)
    tempName = ['temp_' num2str(n) '.fas'];
    filePathList{end+1,1} = tempName;
    fastawrite(tempName, Header((n-1)*seqNum+1:n*seqNum), Sequence((n-1)*seqNum+1:n*seqNum));
end

if floor(length(Header)/seqNum) < ceil(length(Header)/seqNum)
    n = ceil(length(Header)/seqNum);
    tempName = ['temp_' num2str(n) '.fas'];
    filePathList{end+1,1} = tempName;
    fastawrite(tempName, Header((n-1)*seqNum+1:end), Sequence((n-1)*seqNum+1:end));
end

