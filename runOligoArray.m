function runOligoArray(species,varargin)

if length(varargin) >= 1
    params = varargin{1};
else
    params = struct('species','Mouse','verbose',1,'GbMem',16,'maxFragment',1E3,...
        'crosshybeT',75,'secstructT',75.0,'minTm',75,'maxTm',87,...
        'minPercentGC',30.0,'maxPercentGC',50.0,'probeLength',30,...
        'probeLengthMax',30,'minProbeDist',40,...
        'maskedSeq','"GGGGGG;CCCCCC;TTTTTT;AAAAAA"',...
        'oligoArrayPath','C:\OligoArray\','oligoArrayExe','C:\OligoArray\OligoArray2.jar',...
        'blastDbPath','C:\OligoArray\BlastDb\','numParallel',10);
end
params(1).species=species;

%% Set up files for OligoArray; Not decided whether ncrna should be included in the database
if params.verbose
    disp('setting up OligoArray run');
end

if exist([params.oligoArrayPath params.species '.target1kb.fas'], 'file')
    delete([params.oligoArrayPath params.species '.target1kb.fas']);
end
copyfile([params.species '.target1kb.txt'], [params.oligoArrayPath params.species '.target1kb.fas']);

if exist([params.blastDbPath params.species '.cdna.fas'], 'file')
    delete([params.blastDbPath params.species '.cdna.fas*']);
end
copyfile([params.species '.cdna.fas'], [params.blastDbPath params.species '.cdna.fas']);
blastformat('Inputdb', [params.blastDbPath params.species '.cdna.fas'],...
    'FormatArgs', '-o T -p F');

if exist([params.oligoArrayPath 'oligos.txt'], 'file')
    delete([params.oligoArrayPath 'oligos.txt']);
end
if exist([params.oligoArrayPath 'rejected.fas'], 'file')
    delete([params.oligoArrayPath 'rejected.fas']);
end
if exist([params.oligoArrayPath 'OligoArray.log'], 'file')
    delete([params.oligoArrayPath 'OligoArray.log']);
end

%% Execute OligoArray
oligoArrayCommand = ['java -Xmx',ceil(num2str(params.GbMem)),'g -jar ' params.oligoArrayExe ...
    ' -i ' params.oligoArrayPath params.species '.target1kb.fas' ...
    ' -d ' params.blastDbPath params.species '.cdna.fas'...
    ' -o ' params.oligoArrayPath 'oligos.txt' ...
    ' -r ' params.oligoArrayPath 'rejected.fas' ...
    ' -R ' params.oligoArrayPath 'OligoArray.log' ...
    ' -n ' num2str(round(params.maxFragment/params.minProbeDist)) ... % The maximum number of oligos per target sequence
    ' -l ' num2str(params.probeLength) ... % Minimum oligo length
    ' -L ' num2str(params.probeLengthMax) ... % Maximum oligo length
    ' -D ' num2str(params.maxFragment) ... % Maximum distance between 5' of oligo and 3' of target sequence
    ' -t ' num2str(params.minTm) ... % Minimum Tm
    ' -T ' num2str(params.maxTm) ... % Maximum Tm
    ' -s ' num2str(params.secstructT) ... % Threshold temperature for secondary structure prediction
    ' -x ' num2str(params.crosshybeT) ... % Threshold for cross-hybridization
    ' -p ' num2str(params.minPercentGC) ... % Minimum GC content
    ' -P ' num2str(params.maxPercentGC) ... % Maximum GC content
    ' -m ' params.maskedSeq ... % Masked sequences
    ' -g ' num2str(params.minProbeDist) ... % Minimum distance between 5' end of oligos
    ' -N ' num2str(params.numParallel)]; % Number of processes run in parallel

system(oligoArrayCommand);
if params.verbose
    disp('done generating oligos');
end

tempOligos=[params.species '.tempoligos.txt'];
if exist(tempOligos,'file')
    delete(tempOligos);
end
copyfile([params.oligoArrayPath 'oligos.txt'], tempOligos);

