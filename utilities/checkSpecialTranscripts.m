function notSpecialTranscript = checkSpecialTranscripts(transcriptName,params)

% params = struct('species','Mouse','verbose',1,...
%     'number',48,'seqNum',1000,'thres',30,'querySize',30,...
%     'DbSize',2*10^5,'blastArgs','-S 2','parallel', 0,...
%     'specialTranscripts','C:\FISHerMan\Db\Mouse.STList.fas');

[specialTranscripts,~]=fastaread(params(1).specialTranscripts);
specialTranscripts=specialTranscripts';

% transcriptName is one-segment Header
notSpecialTranscript = true;
for n = 1:length(specialTranscripts)
    notSpecialTranscript = notSpecialTranscript && ...
        isempty(strfind(transcriptName,specialTranscripts{n,1}));
end
