function [expressedHeader, expressedSequence] = pickExpressedSeq(seqData, Header, Sequence)

if isempty(seqData{1,1})
    transcriptID = seqData(:,2);

    pos = regexp(Header, ':');
    temp = Header;
    index = [];
    for n = 1:length(Header)
        if length(pos{n,1}) == 2
            temp{n,1} = Header{n,1}(pos{n,1}(1)+1:pos{n,1}(2)-1);
        elseif length(pos{n,1}) == 1
            temp{n,1} = Header{n,1}(pos{n,1}(1)+1:end);
        else
            disp('missing gene ID');
            quit;
        end
        if sum(strcmp(transcriptID, temp{n,1})) ~= 0
            index = [index;n];
        end
    end

    expressedHeader = Header(index);
    expressedSequence = Sequence(index);
else
    transcriptID = seqData(:,1);

    pos = regexp(Header, ':');
    temp = Header;
    for n = 1:length(Header)
        if ~isempty(pos{n,1})
            temp{n,1} = Header{n,1}(1:pos{n,1}(1)-1);
        else
            temp{n,1} = Header{n,1};
        end
    end

    [~,index,~] = intersect(temp, transcriptID, 'stable');
    expressedHeader = Header(index);
    expressedSequence = Sequence(index);
end
