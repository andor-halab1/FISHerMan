%% Before this step, run a portion of main.m to get params
function combineTwoLibs(FISHerManPath,params)

cd([FISHerManPath 'accessories']);
params.onePCR.species = 'combined';
params.probeList.species = 'combined';

lib1 = input('where is the first library you want to combine: ');
lib2 = input('where is the second library you want to combine: ');

[probeHeader1,probeSequence1]=fastaread(lib1);
probeHeader1=probeHeader1';
probeSequence1=probeSequence1';
[probeHeader2,probeSequence2]=fastaread(lib2);
probeHeader2=probeHeader2';
probeSequence2=probeSequence2';

probeHeader = vertcat(probeHeader1,probeHeader2);
probeSequence = vertcat(probeSequence1,probeSequence2);

adapterHeader = probeHeader;
adapterSequence = probeSequence;
probeSequence2Seg = probeSequence;
probeSequenceCore = probeSequence;
for n = 1:length(probeHeader)
    adapterSequence{n,1} = probeSequence{n,1}(21:40);
    probeSequence2Seg{n,1} = probeSequence{n,1}(21:68);
    probeSequenceCore{n,1} = probeSequence{n,1}(31:75);
end

adapterHeader = unique(adapterHeader);
adapterSequence = unique(adapterSequence);

if length(adapterHeader) ~= length(probeHeader)/params.probeList.number
    disp('same transcripts appeared in both libraries');
    quit;
end

adapterList='combined.adapters.txt';
if exist(adapterList, 'file')
    delete(adapterList);
end
fastawrite(adapterList,adapterHeader,adapterSequence);

[probeHeader,probeSequence,probeSequence2Seg,probeSequenceCore]...
    =blast1stPCR(adapterList,probeHeader,probeSequence,probeSequence2Seg,probeSequenceCore,params.onePCR);

params.probeList.number=params.probeList.number-3;
probeList=generateProbeList(adapterList,probeHeader,probeSequence,params.probeList);
