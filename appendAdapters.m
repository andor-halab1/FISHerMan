function [adapterList,probeHeader,probeSequence,probeSequence3Seg,probeSequenceCore]...
    =appendAdapters(adapterList,oligos,params)

% if length(varargin) >= 1
%     params = varargin{1};
% else
%     params = struct('species','Mouse','verbose',1,'keys','ENS\w*T\d*',...
%         'gf','GGAATCGTTGCGGGTGTCCT','grr','CCGCAACATCCAGCATCGTG',...
%         'T7r','CCCTATAGTGAGTCGTATTA',...
%         'rRr','AGAGTGAGTAGTAGTGGAGT','rGr','GATGATGTAGTAGTAAGGGT',...
%         'rBr','TGTGATGGAAGTTAGAGGGT','rIRr','GGAGTAGTTGGTTGTTAGGA');
% end

if params(1).verbose
    disp('concatenating oligos with adapters');
end

[Header, Sequence] = fastaread(oligos);
Header = Header';
Sequence = Sequence';

pos = regexp(Header, params(1).keys, 'end');
trimmedHeader = Header;
for n = 1:length(Header)
    trimmedHeader{n,1} = Header{n,1}(1:pos{n,1});
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

% adapterRSequence = adapterSequence;
% adapterGSequence = adapterSequence;
% adapterBSequence = adapterSequence;
% adapterIRSequence = adapterSequence;
% for n = 1:length(adapterSequence)
%     adapterRSequence{n,1} = strcat(params.rRr, adapterSequence{n,1});
%     adapterGSequence{n,1} = strcat(params.rGr, adapterSequence{n,1});
%     adapterBSequence{n,1} = strcat(params.rBr, adapterSequence{n,1});
%     adapterIRSequence{n,1} = strcat(params.rIRr, adapterSequence{n,1});
% end

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

