function notSpecialTranscript = checkSpecialTranscripts(transcriptName,params)

% params = struct('species','Mouse','verbose',1,...
%     'specialTranscripts','C:\FISHerMan\Db\Mouse.STList.fas');

[specialTranscripts,~]=fastaread(params(1).specialTranscripts);
specialTranscripts=specialTranscripts';

notSpecialTranscript = true;
for n = 1:length(specialTranscripts)
    notSpecialTranscript = notSpecialTranscript && ...
        isempty(strfind(transcriptName,specialTranscripts{n,1}));
end
