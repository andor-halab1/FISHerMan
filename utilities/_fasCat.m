function seq_cat = fasCat(seq, header, footer)

seq_cat = seq;
for n = 2:2:length(seq)
    probe = strcat(header, seq{n,1}, footer);
    seq_cat{n,1} = probe;
end

