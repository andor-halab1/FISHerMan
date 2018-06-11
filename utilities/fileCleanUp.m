function fileCleanUp(filePathList)

    for n = 1:length(filePathList)
        delete(filePathList{n,1});
    end

