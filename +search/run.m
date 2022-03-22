function seq = run(finalbody, launchdates, opts)

%% ARGUMENT VALIDATION
arguments
    finalbody
    launchdates
    opts.FlybyBodies = 2:5
    opts.IterationBudget = 1e5
    opts.NumChildren = 10
    opts.dvBudget = 10
end

opts.FinalBody = finalbody;

%% RUNNING SEARCH
% CREATING TREE ROOT
tree = struct('LaunchNodes', node);
for i = 1:length(launchdates)
    tree.LaunchNodes(i) = node(3, launchdates(i), 0, [], i);
end

for i = 1:1%opts.IterationBudget
    noi = select(tree, opts);
%     if iscell(noi)
%         tree = editProps
%     end
end

%% EXPORTING
seq = tree;

end