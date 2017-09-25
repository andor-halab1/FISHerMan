[~, temp] = xlsread('C:\FISHerMan\Db\refseq2ensembl.xlsx');

data = temp(2:end,:);

for n = 1:length(data)
    temp=regexp(data{n,2},'N\w\w\d*', 'end');
    data{n,2}=data{n,2}(1:temp);
end

missed = [];
for n = 1:length(seqData)
%     temp=regexp(seqData{n,1},'N\w\w\d*', 'end');
%     seqData{n,1}=seqData{n,1}(1:temp);
    index=strcmp(data(:,2),seqData{n,1});
    if sum(index) == 1
        seqData{n,1}=data{index,3};
        seqData{n,2}=data{index,1};
    else
        missed=[missed;n];
    end
end

missedSeq=seqData(missed,:);
seqData(missed,:)=[];

xlswrite('C:\FISHerMan\Db\temp.xlsx',seqData);
