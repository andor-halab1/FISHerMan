probeList = 'C:\FISHerMan\Mouse\Mouse.probes.selected.fas';
[probeHeader,probeSequence]=fastaread(probeList);
probeHeader=probeHeader';
probeSequence=probeSequence';

indexTotal = randperm(length(probeHeader))';
probeHeaderRandom = probeHeader(indexTotal);
probeSequenceRandom = probeSequence(indexTotal);
if exist(probeList, 'file')
    delete(probeList);
end
fastawrite(probeList,probeHeaderRandom,probeSequenceRandom);

% [~,index,indexRandom]=intersect(probeSequence,probeSequenceRandom);
