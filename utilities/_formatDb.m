function formatDb(blastDb)

commandtext=['formatdb.exe -i ' blastDb ' -o T -p F'];
system(commandtext);
