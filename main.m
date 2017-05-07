%% Specify the species to design FISH probes against
disp('weather is good for FISHing! FISHerMan is at work');
species = input('type in the species you would like to design FISH probes against: ');
if isempty(species)
    species = 'Mouse';
end

disp('setting up directories');
FISHerManPath = 'C:\FISHerMan\';
cd([FISHerManPath species]);
addpath(FISHerManPath);
addpath([FISHerManPath 'utilities']);

cdna = input('input the directory where the cdna file can be found: ');
ncrna = input('input the directory where the ncrna file can be found: ');
trna = input('input the directory where a separate trna file can be found: ');
seqData1 = input('input the directory where the RNA-seq data file can be found: ');
if ~isempty(seqData1)
    answer = input('take the average of two RNA-seq replicates? (1/0) ');
    if answer == 1
        seqData2 = input('input the directory where the 2nd RNA-seq data file can be found: ');
        seqData = averageRNASeq(seqData1,seqData2);
    else
        seqData = readRNASeq(seqData1);
    end
else
    seqData = seqData1;
end
transcriptList = input('input the directory where the list of target transcripts can be found: ');
adapterList = input('input the directory where the list of adapters can be found: ');

%% Process the input files
%  now I am just using the mRNA seq data, and it tells me which transcripts
%  in the cdna databse are expressed. Later on, I can also use the total
%  RNA seq data, and it will tell me which transcripts in the ncrna
%  database are expressed. But be sure to include rRNA and tRNA, for often
%  these two types of RNA are depleted in RNA seq.
if isempty(seqData)
    [cdna,cdnaHeader,cdnaSequence]=cdnaParse(cdna);
    [ncrna,trna,ncrnaHeader,ncrnaSequence]=ncrnaParse(ncrna,trna);
else
    [cdna,cdnaHeader,cdnaSequence]=cdnaParse(cdna,seqData);
    [ncrna,trna,ncrnaHeader,ncrnaSequence]=ncrnaParse(ncrna,trna);
end

[abundantrna,abundantrnaHeader,abundantrnaSequence]...
    =abundantrnaParse(cdnaHeader,cdnaSequence,ncrnaHeader,ncrnaSequence,seqData);

[transcriptList,transcriptHeader,transcriptSequence]...
    =transcriptListParse(transcriptList,cdnaHeader,cdnaSequence,ncrnaHeader,ncrnaSequence);

%% Run OligoArray to generate a raw list of oligos
runOligoArray;
oligoList=oligosParse;

%% Append pre-designed adapters to the raw list of oligos
[adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =appendAdapters(adapterList,oligoList);

%% Remove probes that non-specifically bind to primers in the 1st PCR step
[probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =blast1stPCR(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore);

%% Save the probes of each transcripts into individual files
[probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =blastOtherSteps(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore);

%% Remove non-specific probes against the abundant rna database
[probeHeader,probeSequence3Seg,probeSequence,probeSequenceCore]...
    =blastAbundantRNA(adapterList,probeHeader,probeSequence3Seg,probeSequence,probeSequenceCore);

%% Generate the list of probes
probeList=generateProbeList(adapterList,probeHeader,probeSequence);

disp('done designing FISH probes');
disp('FISHerMan is at rest');

