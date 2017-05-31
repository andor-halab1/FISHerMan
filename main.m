%% Specify the species to design FISH probes against
disp('weather is good for FISHing! FISHerMan is at work');

species = input('type in the species you would like to design FISH probes against: ');
if isempty(species)
    species = 'Mouse';
end

disp('setting up directories');
FISHerManPath = 'C:\FISHerMan\';
if exist([FISHerManPath species],'dir')
    cd([FISHerManPath species]);
else
    mkdir([FISHerManPath species]);
    cd([FISHerManPath species]);
end
addpath(FISHerManPath);
addpath([FISHerManPath 'utilities']);

parameters = input('input the directory where the parameter file can be found: ');
params=readParameters(species,parameters);

%% Process the RNA-Seq files
if params.rnaSeq(1).data == 0
    seqData = [];
elseif params.rnaSeq(1).data == 1
    seqData = readRNASeq(params.rnaSeq(1).dir1,params.rnaSeq);
elseif params.rnaSeq(1).data == 2
    seqData = averageRNASeq(params.rnaSeq(1).dir1,params.rnaSeq(1).dir2,params.rnaSeq);
end

%% Process the database files
%  now I am just using the mRNA seq data, and it tells me which transcripts
%  in the cdna databse are expressed. Later on, I can also use the total
%  RNA seq data, and it will tell me which transcripts in the ncrna
%  database are expressed. But be sure to include rRNA and tRNA, for often
%  these two types of RNA are depleted in RNA seq.
[cdnaHeader,cdnaSequence]=cdnaParse(params.cdna(1).dir1,seqData,params.cdna);
if params.rnaSeq(1).mRNA
    [ncrnaHeader,ncrnaSequence]...
        =ncrnaParse(params.ncrna(1).dir1,[],params.ncrna(1).dirT,params.ncrna);
else
    [ncrnaHeader,ncrnaSequence]...
        =ncrnaParse(params.ncrna(1).dir1,seqData,params.ncrna(1).dirT,params.ncrna);
end

[abundantrnaHeader,abundantrnaSequence]...
    =abundantrnaParse(cdnaHeader,cdnaSequence,ncrnaHeader,ncrnaSequence,seqData,params.abundantrna);

%% Process the transcript list file
[transcriptHeader,transcriptSequence]...
    =transcriptListParse(params.transcriptList(1).dir1,cdnaHeader,cdnaSequence,ncrnaHeader,ncrnaSequence,params.transcriptList);

%% Run OligoArray to generate a raw list of oligos
runOligoArray;
oligoList=oligosParse(params.oligos);

%% Append pre-designed adapters to the raw list of oligos
[adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =appendAdapters(params.adapters(1).dir1,oligoList,params.adapters);

%% Remove probes that non-specifically bind to primers in the 1st PCR step
[probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =blast1stPCR(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore,params.onePCR);

%% Remove non-specific probes that will affect other synthesis steps
[probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =blastOtherSteps(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore,params.otherSteps);

%% Remove non-specific probes against the abundant rna database
[probeHeader,probeSequence3Seg,probeSequence,probeSequenceCore]...
    =blastAbundantRNA(adapterList,probeHeader,probeSequence3Seg,probeSequence,probeSequenceCore,params.arna);

%% Generate the list of probes
probeList=generateProbeList(adapterList,probeHeader,probeSequence,params.probeList);

disp('done designing FISH probes');
disp('FISHerMan is at rest');

