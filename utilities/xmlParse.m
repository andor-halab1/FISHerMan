function branch = xmlParse(root, rootName, branchName)

    try
        branch = root.getElementsByTagName(branchName);
        branch = branch.item(0);
    catch
        disp(['error occurred at ' rootName]);
    end
