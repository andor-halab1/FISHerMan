function seqData = averageRNASeq(seqData1,seqData2,params)

% params = struct('species','Mouse','verbose',1,'data',2,...
%     'dir1','C:\FISHerMan\Db\mouse_frontal_cortex_mRNASeq_ENCFF653BKJ.xlsx',...
%     'dir2','C:\FISHerMan\Db\mouse_frontal_cortex_mRNASeq_ENCFF703SOK.xlsx',...
%     'mRNA',1,'keys',{'ENS\w*T\d*','ENS\w*G\d*'},...
%     'thres',0.1);

seqData1 = readRNASeq(seqData1,params);
seqData2 = readRNASeq(seqData2,params);

transcriptID1 = seqData1(:,1);
transcriptID2 = seqData2(:,1);

[~,index1,index2] = intersect(transcriptID1, transcriptID2, 'stable');

seqData1 = seqData1(index1,:);
seqData2 = seqData2(index2,:);

seqData = seqData1;
for n = 1:length(seqData)
    if strcmp(seqData1{n,1},seqData2{n,1})
        seqData{n,3} = (seqData1{n,3}+seqData2{n,3})/2;
    else
        warning('transcript ID not matched');
    end
end

FPKM = cell2mat(seqData(:,3));
index = FPKM>=params(1).thres;
seqData = seqData(index,:);
