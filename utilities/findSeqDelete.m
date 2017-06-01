function seqDelete=findSeqDelete(data)

%% Delete probes that non-specifically bind 2nd PCR primers
seqDelete = [];
for n = 1:length(data{2,1})
    if ~isempty(data{2,1}(n).Hits)
        seqDelete = [seqDelete n];
    end
end

%% Delete probes that non-specifically bind other probes
pairDelete = {};
for n = 1:length(data{1,1})
    for m = 1:length(data{1,1}(n).Hits)
        pairDelete{end+1,1} = data{1,1}(n).Query;
        pairDelete{end,2} = data{1,1}(n).Hits(m).Name;
    end
end

pos = regexp(pairDelete,'=');
pairNumber = {};
for k = 1:size(pos,1)
    temp1 = pairDelete{k,1}(pos{k,1}+1:end);
    temp2 = pairDelete{k,2}(pos{k,2}+1:end);
    temp1 = uint16(str2double(temp1));
    temp2 = uint16(str2double(temp2));
    pairNumber{end+1,1} = temp1;
    pairNumber{end,2} = temp2;
end
pairNumber = cell2mat(pairNumber);

while size(pairNumber,1) > 0
    temp = reshape(pairNumber, [numel(pairNumber),1]);
    tempDelete = mode(temp);
    seqDelete = [seqDelete tempDelete];
    index = (pairNumber == tempDelete);
    index = logical(index(:,1)+index(:,2));
    pairNumber(index,:) = [];
end

seqDelete = unique(seqDelete, 'stable');
seqDelete = seqDelete';

