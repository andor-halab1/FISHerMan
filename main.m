%% Specify the species to design FISH probes against
disp('a good weather for FISHing! FISHerMan is at work');
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
seqData = input('input the directory where the RNA-seq data file can be found: ');
if ~isempty(seqData)
    answer = input('take the average of two RNA-seq replicates? (1/0) ');
    if answer == 1
        seqData2 = input('input the directory where the 2nd RNA-seq data file can be found: ');
        seqData = averageSeqData(seqData,seqData2);
    end
end
transcriptList = input('input the directory where the list of target transcripts can be found: ');

%% Process the input files
if isempty(seqData)
    [cdna,cdnaHeader,cdnaSequence]=cdnaParse(cdna);
    [ncrna,ncrnaHeader,ncrnaSequence]=ncrnaParse(ncrna);
else
    if answer ~= 1
        seqData = readRNASeq(seqData);
    end
    [cdna,cdnaHeader,cdnaSequence]=cdnaParse(cdna,seqData);
    [ncrna,ncrnaHeader,ncrnaSequence]=ncrnaParse(ncrna);
end

[abundantrna,abundantrnaHeader,abundantrnaSequence]...
    =abundantrnaParse(cdnaHeader,cdnaSequence,ncrnaHeader,ncrnaSequence,seqData);

[transcriptList,transcriptHeader,transcriptSequence]...
    =transcriptListParse(transcriptList,cdnaHeader,cdnaSequence,ncrnaHeader,ncrnaSequence);

%% Run OligoArray to generate a raw list of oligos
oligoList = runOligoArray;

%% Append pre-designed adapters to the raw list of oligos
adapterList = input('input the directory where the list of adapters can be found: ');
[adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =adapterAppend(adapterList,oligoList);

%% Remove probes that non-specifically bind to primers in the 1st PCR step
[probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =blast1stPCR(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore);

%% Save the probes of each transcripts into individual files
[probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =blastOtherSteps(adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore);

%% Remove non-specific probes against the abundant rna database
[probeHeader,probeSequence3Seg,probeSequence,probeSequenceCore]...
    =blastAbundantRNA(adapterList,probeHeader,probeSequence3Seg,probeSequence,probeSequenceCore);

%% Save the list of probes
disp('saving the list of probes');
probeList = [species '.probes.fas'];
if exist(probeList, 'file')
    delete(probeList);
end
fastawrite(probeList,probeHeader,probeSequence);

disp('done designing FISH probes');
disp('FISHerMan is at rest');

