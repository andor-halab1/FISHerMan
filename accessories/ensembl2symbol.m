[longHeader,~]=fastaread('C:\FISHerMan\Db\Mouse.cdna.fa');
longHeader=longHeader';

[pos1{1,1}, pos2{1,1}] = regexp(longHeader, 'ENS\w*T\d*', 'start', 'end');
[pos1{2,1}, pos2{2,1}] = regexp(longHeader, 'gene_symbol:', 'start', 'end');
[pos1{3,1}, pos2{3,1}] = regexp(longHeader, 'description:', 'start', 'end');

for n = 1:length(longHeader)
    entry{n,1}=longHeader{n,1}(1:pos2{1,1}{n,1});
    entry{n,2}=longHeader{n,1}(pos2{2,1}{n,1}+1:pos1{3,1}{n,1}-2);
end

[Header,~]=fastaread('C:\FISHerMan\designed.libraries\combined\combined.adapters.txt');
Header=Header';

for n = 1:length(Header)
    temp=regexp(Header{n,1},':');
    if ~isempty(temp)
        Header{n,1}=Header{n,1}(1:temp(1)-1);
    end
end

[~,~,index]=intersect(Header,entry(:,1),'stable');

matchedHeader=entry(index,1);
matchedSymbol=entry(index,2);

fastawrite('C:\FISHerMan\designed.libraries\combined\combined.symbols.txt',matchedHeader,matchedSymbol);
