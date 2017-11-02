[Header,Sequence]=fastaread('C:\FISHerMan\designed.libraries\MN\MN.probes.nr.txt');
Header=Header';
Sequence=Sequence';

[Header2,Sequence2]=fastaread('C:\FISHerMan\designed.libraries\MN.probes.nr.txt');
Header2=Header2';
Sequence2=Sequence2';

[~,probeSequence]=fastaread('C:\FISHerMan\designed.libraries\MN\MN.oligos.txt');
probeSequence=probeSequence';

[~,probeSequence2]=fastaread('C:\FISHerMan\designed.libraries\MN.oligos.txt');
probeSequence2=probeSequence2';

uniqueHeader=unique(Header,'stable');
uniqueHeader2=unique(Header2,'stable');

[commonHeader,index,index2]=intersect(uniqueHeader,uniqueHeader2,'stable');
[commonSequence,indexS,indexS2]=intersect(probeSequence,probeSequence2,'stable');
