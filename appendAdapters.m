function [adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =appendAdapters(adapterList,oligos,params)

% params = struct('species','Mouse','verbose',1,...
%     'gf','GGAATCGTTGCGGGTGTCCT','grr','CCGCAACATCCAGCATCGTG',...
%     'T7r','CCCTATAGTGAGTCGTATTA',...
%     'rRr','AGAGTGAGTAGTAGTGGAGT','rGr','GATGATGTAGTAGTAAGGGT',...
%     'rBr','TGTGATGGAAGTTAGAGGGT','rIRr','GGAGTAGTTGGTTGTTAGGA');

if params(1).verbose
    disp('concatenating oligos with adapters');
end

[Header, Sequence] = fastaread(oligos);
Header = Header';
Sequence = Sequence';

pos = regexp(Header, ':');
trimmedHeader = Header;
for n = 1:length(Header)
    trimmedHeader{n,1} = Header{n,1}(1:pos{n,1}(1)-1);
end
uniqueHeader = unique(trimmedHeader, 'stable');

geneNum = length(uniqueHeader);
[~, adapterSequence] = fastaread(adapterList);
adapterHeader = uniqueHeader;
adapterSequence = adapterSequence(1:geneNum)';

adapterList = [params(1).species '.adapters.fas'];
if exist(adapterList, 'file')
    delete(adapterList);
end
fastawrite(adapterList, adapterHeader, adapterSequence);

probeHeader = {};
probeSequence = {};
probeSequence3Seg = {};
probeSequenceCore = {};
for n = 1:length(adapterHeader)
    if params(1).verbose && mod(n, 1000) == 1
        disp(['  concatenating oligos for transcript no. ' num2str(n)]);
    end
    index = ismember(trimmedHeader, adapterHeader{n,1});
    probeHeader(end+1:end+sum(index),1) = Header(index);
    probe = Sequence(index);
    probe3Seg = Sequence(index);
    probeCore = Sequence(index);
    for m = 1:length(probe)
        temp = probe{m,1};
        temp3Seg = probe3Seg{m,1};
        tempCore = probeCore{m,1};
        
        temp = strcat(params(1).gf, adapterSequence{n,1}, temp, params(1).grr);
        temp3Seg = strcat(adapterSequence{n,1}, temp3Seg, params(1).grr, 'CCC');
        tempCore = strcat(adapterSequence{n,1}(11:20), tempCore, params(1).grr(1:10));
        
        probe{m,1} = temp;
        probe3Seg{m,1} = temp3Seg;
        probeCore{m,1} = tempCore;
    end
    probeSequence(end+1:end+sum(index),1) = probe;
    probeSequence3Seg(end+1:end+sum(index),1) = probe3Seg;
    probeSequenceCore(end+1:end+sum(index),1) = probeCore;
end

% if params(1).verbose
%     disp('saving the probe fasta file');
% end
% 
% probes = [params.species '.probes.fas'];
% if exist(probes, 'file')
%     delete(probes);
% end
% fastawrite(probes, probeHeader, probeSequence);

