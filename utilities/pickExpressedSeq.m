function [expressedHeader, expressedSequence] = pickExpressedSeq(seqData, Header, Sequence)

transcriptID = seqData(:,1);

pos = regexp(Header, ':');
temp = Header;
for n = 1:length(Header)
    temp{n,1} = Header{n,1}(1:pos{n,1}(1)-1);
end

[~,index,~] = intersect(temp, transcriptID, 'stable');
expressedHeader = Header(index);
expressedSequence = Sequence(index);
