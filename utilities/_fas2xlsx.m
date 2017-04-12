[Header, Sequence] = fastaread('C:\Users\Rong Li Lab\Documents\MATLAB\bc25mer.240k.fasta');

Header = Header';
Sequence = Sequence';
xlswrite('C:\Users\Rong Li Lab\Documents\MATLAB\bc25mer.240k.xlsx',Sequence);

% fastawrite('C:\Users\Rong Li Lab\Documents\MATLAB\bc25mer.240k.fas', Header, Sequence);


