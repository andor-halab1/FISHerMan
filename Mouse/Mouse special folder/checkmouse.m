%% check if gene-specific primers are correct

[Header,Sequence]=fastaread('C:\FISHerMan\Mouse\Mouse special folder\Mouse.probestemp.txt');
Header=Header';
Sequence=Sequence';

tempHeader=Header;
tempSequence=Sequence;
for n = 1:length(Sequence)
    tempSequence{n,1}=Sequence{n,1}(21:40);
end

[totalHeader,totalSequence]=fastaread('C:\FISHerMan\Db\Mouse.alladapters.txt');
totalHeader=totalHeader';
totalSequence=totalSequence';

uniqueTempSequence=unique(tempSequence,'stable');

[~,index1,index2]=intersect(totalSequence,uniqueTempSequence);
index2=sort(index2);

commonTempSequence=uniqueTempSequence(index2);

indexTotal=ismember(tempSequence,commonTempSequence{1,1});
for n = 2:length(commonTempSequence)
    index=ismember(tempSequence,commonTempSequence{n,1});
    indexTotal=indexTotal+index;
end
indexTotal=logical(indexTotal);

Header=Header(indexTotal);
Sequence=Sequence(indexTotal);

fastawrite('C:\FISHerMan\Mouse\Mouse special folder\Mouse.probes.temp.txt',Header,Sequence);

%% check non-specific binding against other mRNA

[Header,Sequence]=fastaread('C:\FISHerMan\Mouse\Mouse special folder\Mouse.probes.temp.txt');
Header=Header';
Sequence=Sequence';

tempHeader=Header;
tempSequence=Sequence;
for n = 1:length(Sequence)
    tempSequence{n,1}=Sequence{n,1}(41:70);
end

[totalHeader,totalSequence]=fastaread('C:\FISHerMan\Mouse\Mouse special folder\Mouse.oligos.temp.txt');
totalHeader=totalHeader';
totalSequence=totalSequence';

[~,index1,index2]=intersect(totalSequence,tempSequence);
index2=sort(index2);

Header=Header(index2);
Sequence=Sequence(index2);

uniqueHeader=unique(Header);

indexTotal=ismember(Header,uniqueHeader{1,1});
for n = 2:length(uniqueHeader)
    index=ismember(Header,uniqueHeader{n,1});
    if sum(index)<48 && isempty(strfind(uniqueHeader{n,1},'ENSSPT'))
        indexTotal=indexTotal+index;
    end
end
indexTotal=logical(indexTotal);

Header(indexTotal)=[];
Sequence(indexTotal)=[];

fastawrite('C:\FISHerMan\Mouse\Mouse special folder\Mouse.probes.mRNA.txt',Header,Sequence);

%% start over to blast against mouse abundant rna

[Header,Sequence]=fastaread('C:\FISHerMan\Mouse\Mouse special folder\Mouse.probes.mRNA.txt');
Header=Header';
Sequence=Sequence';

tempHeader=Header;
tempSequence=Sequence;
for n = 1:length(Sequence)
    tempSequence{n,1}=strcat(Sequence{n,1}(21:90),'CCC');
end

[Header,tempSequence,Sequence]...
    =blastAbundantRNASimple(Header,tempSequence,Sequence,params.arna);

uniqueHeader=unique(Header);

% this number (10) was tested one-by-one from 1
indexTotal=ismember(Header,uniqueHeader{10,1});
for n = 2:length(uniqueHeader)
    index=ismember(Header,uniqueHeader{n,1});
    if sum(index)<48 && isempty(strfind(uniqueHeader{n,1},'ENSSPT'))
        indexTotal=indexTotal+index;
    end
end
indexTotal=logical(indexTotal);

Header(indexTotal)=[];
Sequence(indexTotal)=[];

fastawrite('C:\FISHerMan\Mouse\Mouse special folder\Mouse.probes.arna.txt',Header,Sequence);




