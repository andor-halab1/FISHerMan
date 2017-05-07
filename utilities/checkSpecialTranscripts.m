function notSpecialTranscript = checkSpecialTranscripts(transcriptName,varargin)

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,...
        'specialTranscripts',{'ENSSPT','ENSMUST00000100497','ENSMUST00000118875'});
end

notSpecialTranscript = true;
for n = 1:length(params)
    notSpecialTranscript = notSpecialTranscript && ...
        isempty(strfind(transcriptName,params(n).specialTranscripts));
end
