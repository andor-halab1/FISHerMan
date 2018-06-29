function data = blastOp(filePath, DbPath, blastArgs)
    data = blastlocal('InputQuery', filePath, 'Database', DbPath, ...
        'Program', 'blastn', 'BLASTArgs', blastArgs);
end
