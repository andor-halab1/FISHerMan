[dum, temp] = xlsread('C:\Users\Rong Li Lab\Downloads\GSE62571_RAW\mouse_hippocampus_mRNASeq_1.xlsx');

seqData = cell(size(dum,1),3);
seqData(:,1:2) = temp(2:end,1:2);

for n = 1:size(dum,1)
    pos_t = regexp(seqData{n,1}, 'N\w\w\d*', 'end');
    seqData{n,1} = seqData{n,1}(1:pos_t);
    seqData{n,3} = dum(n,6);
end

FPKM = cell2mat(seqData(:,3));
index = FPKM>=0.4;
seqData = seqData(index,:);

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

[~,index,~]=unique(seqData(:,1));
seqDataUnique=seqData(index,:);

xlswrite('C:\FISHerMan\Db\mouse_hippocampus_mRNASeq_01.xlsx',seqDataUnique);
xlswrite('C:\FISHerMan\Db\mouse_hippocampus_mRNASeq_01_missed.xlsx',missedSeq);
