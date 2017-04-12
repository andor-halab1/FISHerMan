function testDbSize = getDbSize(testDb,varargin)

% testDb = 'C:\FISHerMan\Mouse\Mouse.abundantrnaDb.fas';

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,...
        'thres',40,'querySize',20,'DbSize',5*10^7,'blastArgs','-S 1');
end

eValue = bitScore2eValue(params.thres, params.querySize, params.DbSize);

DbPath = testDb;
blastArgs = [params.blastArgs ' -e ' num2str(eValue)];

[Header, Sequence] = fastaread(testDb);
Header = Header';
Sequence = Sequence';

tempHeader{1,1} = Header{1,1};
tempSequence{1,1} = Sequence{1,1}(1:20);
fastawrite('testDb.fas',tempHeader,tempSequence);

data = blastOp('testDb.fas', DbPath, blastArgs);

testDbSize = 2*data(1).Hits(1).HSPs(1).Expect/(params.querySize*2^(-data(1).Hits(1).HSPs(1).Score));

fileCleanUp({'testDb.fas'});

