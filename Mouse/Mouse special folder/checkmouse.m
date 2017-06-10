[Header,Sequence]=fastaread('C:\FISHerMan\Mouse\Mouse special folder\Mouse.probes.txt');
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
index2=unique(index2);

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

fastawrite('Mouse.probes.temp.txt',Header,Sequence);

%% start over to blast against mouse abundant rna
[Header,Sequence]=fastaread('C:\FISHerMan\Mouse\Mouse special folder\Mouse.probes.temp.txt');
Header=Header';
Sequence=Sequence';

tempHeader=Header;
tempSequence=Sequence;
for n = 1:length(Sequence)
    tempSequence{n,1}=Sequence{n,1}(21:90);
end

[Header,tempSequence,Sequence]...
    =blastAbundantRNASimple(Header,tempSequence,Sequence,params.arna);

uniqueHeader=unique(Header);

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




