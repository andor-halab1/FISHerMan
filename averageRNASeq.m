function seqData = averageRNASeq(seqData1,seqData2,varargin)

% seqData1 = 'mouse_frontal_cortex_mRNASeq_ENCFF653BKJ.xlsx';
% seqData2 = 'mouse_frontal_cortex_mRNASeq_ENCFF703SOK.xlsx';

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,...
        'keys',{'ENS\w*T\d*','ENS\w*G\d*'},'thres',0.1);
end

seqData1 = readRNASeq(seqData1,params);
seqData2 = readRNASeq(seqData2,params);

transcriptID1 = seqData1(:,1);
transcriptID2 = seqData2(:,1);

[~,index1,index2] = intersect(transcriptID1, transcriptID2, 'stable');

seqData1 = seqData1(index1,:);
seqData2 = seqData2(index2,:);

seqData = seqData1;
for n = 1:length(seqData)
    seqData{n,3} = (seqData1{n,3}+seqData2{n,3})/2;
end

FPKM = cell2mat(seqData(:,3));
index = FPKM>=params(1).thres;
seqData = seqData(index,:);
