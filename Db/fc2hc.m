%% Before this step, run a portion of main.m to get seqData
[~,Sequence]=fastaread('C:\FISHerMan\Db\MN.transcriptList.fc.txt');
Sequence=Sequence';

nonExpressed=[];
for n = 1:length(Sequence)
    index=strcmp(seqData(:,1),Sequence{n,1});
    if sum(index) ~= 0
        geneID = seqData{index,2};
        indexGeneID = strcmp(seqData(:,2),geneID);
        rnaSeq = cell2mat(seqData(indexGeneID,3));
        transcriptID = seqData(indexGeneID,1);
        
        [~,temp]=max(rnaSeq);
        Sequence{n,1}=transcriptID{temp,1};
    else
        nonExpressed=[nonExpressed;n];
    end
end

Sequence(nonExpressed,:)=[];

fastawrite('C:\FISHerMan\Db\MN.transcriptList.hc.txt',Sequence,Sequence);
