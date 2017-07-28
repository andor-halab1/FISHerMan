[Header,Sequence]=fastaread('C:\FISHerMan\Mouse\Mouse special folder\Mouse.probestemp.txt');
Header=Header';
Sequence=Sequence';

tempHeader=Header;
tempSequence=Sequence;
for n = 1:length(Sequence)
    tempSequence{n,1}=strcat(Sequence{n,1}(21:90),'CCC');
end

[Header,tempSequence,Sequence]...
    =blastAbundantRNASimple(Header,tempSequence,Sequence,params.arna);

geneDelete = {};
for n = 1:length(data)
    flag = 0;
    for m = 1:length(data(n).Hits)
        if isempty(strfind(data(n).Query,data(n).Hits(m).Name))
            flag = 1;
        end
    end
    if flag == 1
        geneDelete{end+1,1} = data(n).Query;
    end
end

pos=regexp(geneDelete,'=');

for n = 1:length(geneDelete)
    geneDelete{n,1}=geneDelete{n,1}(1:pos{n,1}-1);
end
uniquegeneDelete=unique(geneDelete);

fastawrite('Mouse.probestemp.arna.txt',uniquegeneDelete,uniquegeneDelete);
