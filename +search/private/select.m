function noi = select(tree, inputOpts)

% FINDING LAUNCH NODE OF INTEREST
J = [tree.LaunchNodes(:).J];
[~, idx] = min(J);
noi = tree.LaunchNodes(idx);

% PATHING DOWN TREE
while isfield(noi, 'Children')
    % COLLECTING DATA
    J = [noi.Children(:).J];
    fbd = [noi.Children(:).ID];
    budget = [noi.Children(:).Budget];
    
    % REMOVING COMPLETED NODES
    tf1 = fbd ~= inputOpts.FinalBody;
    
    % REMOVING BUDGET EXCEEDED NODES
    tf2 = budget < inputOpts.dvBudget;
    
    % CULLING VALUES
    idx1 = find(tf1 & tf2);
    
    % SKIPPING IF NO VALID CHILDREN
    if isempty(idx1)
        noi = {NaN, noi.ParentPath};
        return
    end
    
    % CHOOSING LOWEST VALUE
    [~, idx2] = min(J(idx1));
    noi = noi.Children(idx1(idx2));
end
end