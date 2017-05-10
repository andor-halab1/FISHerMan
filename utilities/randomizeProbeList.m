[probeHeader,probeSequence]=fastaread('C:\FISHerMan\Mouse\Mouse.probestemp.fas');
probeHeader=probeHeader';
probeSequence=probeSequence';

indexTotal = randperm(length(probeHeader))';
probeHeaderRandom = probeHeader(indexTotal);
probeSequenceRandom = probeSequence(indexTotal);
probeList = 'C:\FISHerMan\Mouse\Mouse.probes.selected.fas';
if exist(probeList, 'file')
    delete(probeList);
end
fastawrite(probeList,probeHeaderRandom,probeSequenceRandom);

% [~,index,indexRandom]=intersect(probeSequence,probeSequenceRandom);
