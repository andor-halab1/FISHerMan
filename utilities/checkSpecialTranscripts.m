function notSpecialTranscript = checkSpecialTranscripts (transcriptName,params)
    % params = struct('species','Mouse','verbose',1,...
    %     'specialTranscripts','C:\FISHerMan\Db\Mouse.STList.fas');

    [specialTranscripts,~] = fastaread(params(1).specialTranscripts);
    specialTranscripts=specialTranscripts';

    % transcriptName and specialTranscripts are one-segment Header
    notSpecialTranscript = true;
    for n = 1:length(specialTranscripts)
        if strcmp(transcriptName,specialTranscripts{n,1})
            notSpecialTranscript = false
            return
        end
    end
end
