function seqDelete = blastOneTranscript(OTAdapterHeader,OTAdapterSequence,probeHeader,probeSequenceCore,params)
    
    % params = struct('species','Mouse','verbose',1,'keys',...
    %     'thres',22,'querySize',50,...
    %     'blastArgs1','-S 2','blastArgs2','-S 3',...
    %     'grr','CCGCAACATCCAGCATCGTG','T7r','CCCTATAGTGAGTCGTATTA',...
    %     'rRr','AGAGTGAGTAGTAGTGGAGT','rGr','GATGATGTAGTAGTAAGGGT',...
    %     'rBr','TGTGATGGAAGTTAGAGGGT','rIRr','GGAGTAGTTGGTTGTTAGGA');

    % probeHeader and OTHeader are one-segment Header
    index = ismember(probeHeader, OTAdapterHeader);
    temp = 1:length(probeHeader);
    index = temp(index);
    OTHeader = probeHeader(index);
    OTSequenceCore = probeSequenceCore(index);

    %% Generate database files for Blast
    OTDb1 = [params(1).species '.' OTAdapterHeader '.Db1.fas'];
    OTDb2 = [params(1).species '.' OTAdapterHeader '.Db2.fas'];
    % MatLab's use of blastlocal requires short entry names

    for n = 1:length(OTHeader)
        OTHeader{n,1} = strcat(OTHeader{n,1},'=',num2str(n));
    end

    Header{1,1} = 'ENSPRIMERT11';
    Sequence{1,1} = strcat(params(1).rRr,OTAdapterSequence);
    Header{2,1} = 'ENSPRIMERT12';
    Sequence{2,1} = strcat(params(1).rGr,OTAdapterSequence);
    Header{3,1} = 'ENSPRIMERT13';
    Sequence{3,1} = strcat(params(1).rBr,OTAdapterSequence);
    Header{4,1} = 'ENSPRIMERT14';
    Sequence{4,1} = strcat(params(1).rIRr,OTAdapterSequence);
    Header{5,1} = 'ENSPRIMERT15';
    Sequence{5,1} = strcat(params(1).grr,params(1).T7r);

    if exist(OTDb1, 'file')
        delete([OTDb1 '*']);
    end
    fastawrite(OTDb1, OTHeader, OTSequenceCore);
    blastformat('Inputdb', OTDb1,...
        'FormatArgs', '-o T -p F');

    if exist(OTDb2, 'file')
        delete([OTDb2 '*']);
    end
    fastawrite(OTDb2, Header, Sequence);
    blastformat('Inputdb', OTDb2,...
        'FormatArgs', '-o T -p F');

    %% Blast probes against 2nd PCR primers and other probes
    DbPath1 = OTDb1;
    DbPath2 = OTDb2;
    params(1).DbSize1 = getDbSize(DbPath1);
    params(1).DbSize2 = getDbSize(DbPath2);

    eValue1 = bitScore2eValue(params(1).thres, params(1).querySize, params(1).DbSize1);
    eValue2 = bitScore2eValue(params(1).thres, params(1).querySize, params(1).DbSize2);
    blastArgs1 = [params(1).blastArgs1 ' -e ' num2str(eValue1) ' -b ' num2str(length(OTHeader))];
    blastArgs2 = [params(1).blastArgs2 ' -e ' num2str(eValue2) ' -b ' num2str(length(OTHeader))];

    if params(1).verbose
        disp('  blasting the temporary probe list file');
        startTime = tic;
    end
    
    blastData{1,1} = blastOp(OTDb1, DbPath1, blastArgs1);
    blastData{2,1} = blastOp(OTDb1, DbPath2, blastArgs2);
    
    if params(1).verbose
        totalTime = toc(startTime);
        disp(['  elapsed time is ' num2str(totalTime) ' seconds']);
    end

    %% Find out the probes that non-specifically bind 2nd PCR primers and other probes
    seqDelete = findSeqDelete(blastData);

    if ~isempty(seqDelete)
        seqDelete = index(seqDelete);
    end

    delete([OTDb1 '*']);
    delete([OTDb2 '*']);
end
