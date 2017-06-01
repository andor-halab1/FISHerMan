function data = readRNASeq(seqData, params)

% params = struct('species','Mouse','verbose',1,'data',2,...
%     'dir1','C:\FISHerMan\Db\mouse_frontal_cortex_mRNASeq_ENCFF653BKJ.xlsx',...
%     'dir2','C:\FISHerMan\Db\mouse_frontal_cortex_mRNASeq_ENCFF703SOK.xlsx',...
%     'mRNA',1,'keys',{'ENS\w*T\d*','ENS\w*G\d*'},...
%     'thres',0.1);

if params(1).verbose
    disp('reading the RNA-seq data file');
end

[dum, temp] = xlsread(seqData);

data = cell(size(dum,1),3);
data(:,1:2) = temp(2:end,1:2);

for n = 1:size(dum,1)
    pos_t = regexp(data{n,1}, params(1).keys, 'end');
    pos_g = regexp(data{n,2}, params(2).keys, 'end');
    data{n,1} = data{n,1}(1:pos_t);
    data{n,2} = data{n,2}(1:pos_g);
    data{n,3} = dum(n,5);
end

