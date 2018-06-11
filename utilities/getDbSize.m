function testDbSize = getDbSize(testDb,varargin)

    % testDb = 'C:\FISHerMan\Mouse\Mouse.abundantrnaDb.fas';

    if length(varargin) >= 1
        params = varargin{1};
    else
        params = struct('thres', 40, 'querySize', 20, 'DbSize', 1e8, 'blastArgs', '-S 1');
    end

    [Header, Sequence] = fastaread(testDb);
    Header = Header';
    Sequence = Sequence';

    for n = 1:length(Header)
        if length(Sequence{n,1}) >= 20
            break;
        end
    end

    tempHeader{1,1} = Header{n,1};
    tempSequence{1,1} = Sequence{n,1}(1:20);
    fastawrite('testDb.fas',tempHeader,tempSequence);

    eValue = bitScore2eValue(params.thres, params.querySize, params.DbSize);
    DbPath = testDb;
    blastArgs = [params.blastArgs ' -e ' num2str(eValue)];
    data = blastOp('testDb.fas', DbPath, blastArgs);

    if ~isempty(data(1).Hits)
        testDbSize = 2*data(1).Hits(1).HSPs(1).Expect/(params.querySize*2^(-data(1).Hits(1).HSPs(1).Score));
    else
        disp('using a wrong database size');
        testDbSize = 10000;
    end

    fileCleanUp({'testDb.fas'});

