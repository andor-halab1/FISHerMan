function [longHeader,longSequence,X]...
    =transcriptListParse(transcriptList,cdnaHeader,cdnaSequence,ncrnaHeader,ncrnaSequence,params)

% params = struct('species','Mouse','verbose',1,...
%     'dir1','C:\FISHerMan\Db\Mouse.transcriptList.fas',...
%     'length',40,'number',24);

% Pick transcript sequences that are longer than a certain threshold
if params(1).verbose
    disp('reading the name file of target transcripts');
    disp('  thresholding based on sequence length');
end

% use provided transcriptList as targets for designing FISH probes
% targetHeader are one-segment Header
[targetHeader,~]=fastaread(transcriptList);
targetHeader = targetHeader';

totalHeader = vertcat(cdnaHeader,ncrnaHeader);
totalSequence = vertcat(cdnaSequence,ncrnaSequence);
[totalHeader, totalSequence] = pickExpressedSeq(targetHeader, totalHeader, totalSequence);

longHeader = {};
longSequence = {};
for n = 1:length(totalSequence)
    if length(totalSequence{n,1}) >= (params(1).length*params(1).number)
        % this is to make sure longHeader are two-segment Header
        pos=regexp(totalHeader{n,1},':');
        if length(pos) >= 2
            totalHeader{n,1}=totalHeader{n,1}(1:pos(2)-1);
        end
        longHeader{end+1,1} = totalHeader{n,1};
        longSequence{end+1,1} = totalSequence{n,1};
    end
end

if params(1).verbose
    disp(['  ' num2str(length(longHeader)) ' transcripts have length longer than '...
        num2str(params(1).length*params(1).number) ' nts']);
end

%% Check transcript GC contents

GC=[];
for n = 1:length(longHeader)
    G=length(strfind(longSequence{n,1},'G'));
    C=length(strfind(longSequence{n,1},'C'));
    GC=[GC;(G+C)/length(longSequence{n,1})*100];
end

handle = figure(1);
hist(GC,0:5:100);
title('histogram of GC contents');
temp=axis;
temp(1)=0;
temp(2)=100;
axis(temp);

disp('  click to select the boundaries for GC contents');
[X,~]=ginput(2);
close(handle);
X=sort(X,'ascend');
X=floor(X);

if params(1).verbose
    disp(['  GC content threshold is from ' num2str(X(1)) '% to ' num2str(X(2)) '%']);
end

% Segment each transcript sequence into 1kb pieces, so later analysis runs more efficiently
if params(1).verbose
    disp('  segmenting long sequences into 1kb segments');
end

segHeader = {};
segSequence = {};
for n = 1:length(longHeader)
    segN = floor(length(longSequence{n,1})/1000);
    probe = longSequence{n,1};
    for m = 1:segN
        segHeader{end+1,1} = longHeader{n,1};
        segSequence{end+1,1} = probe(1000*(m-1)+1:1000*m);
    end
    if length(probe) >= (1000*segN+30)
        segHeader{end+1,1} = longHeader{n,1};
        segSequence{end+1,1} = probe(1000*segN+1:end);
    end
end

if params(1).verbose
    disp('saving the target transcript fasta file');
end

transcriptList = [params(1).species '.target1kb.txt'];
if exist(transcriptList, 'file')
    delete(transcriptList);
end
fastawrite(transcriptList, segHeader, segSequence);
