[dum, temp] = xlsread('C:\FISHerMan\Db\mouse_frontal_cortex_mRNASeq_ENCFF653BKJ.xlsx');

data = cell(size(dum,1),3);
data(:,1:2) = temp(2:end,1:2);

for n = 1:size(dum,1)
    pos_t = regexp(data{n,1}, 'ENS\w*T\d*', 'end');
    pos_g = regexp(data{n,2}, 'ENS\w*G\d*', 'end');
    data{n,1} = data{n,1}(1:pos_t);
    data{n,2} = data{n,2}(1:pos_g);
    data{n,3} = dum(n,5);
end

geneID=unique(data(:,2),'stable');

for n = 1:length(geneID)
    index=strcmp(data(:,2),geneID{n,1});
    geneID{n,2}=sum(cell2mat(data(index,3)));
end

index = cell2mat(geneID(:,2))<0.5;
geneID(index,:)=[];

[dum, temp] = xlsread('C:\FISHerMan\Db\mouse_hippocampus_totalRNASeq.xlsx');

data2 = cell(size(dum,1),2);
data2(:,1) = temp(2:end,2);
for n = 1:size(dum,1)
    data2{n,2} = dum(n,1);
end

[~,index1,index2]=intersect(geneID(:,1),data2(:,1));

for n = 1:length(index1)
    rnaExpName{n,1}=geneID{index1(n),1};
    rnaExp(n,1)=geneID{index1(n),2};
    rnaExp(n,2)=data2{index2(n),2};
end

figure(1);
loglog(rnaExp(:,1),rnaExp(:,2),'o');
