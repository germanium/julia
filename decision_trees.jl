# From http://bensadeghi.com/decision-trees-julia/
# Random forest example

Pkg.add("DecisionTree");    using DecisionTree
Pkg.add("RDatasets");       using RDatasets

                                        # Load Data
iris = data("datasets", "iris");
features = matrix(iris[:, 2:5]);
labels = vector(iris[:, "Species"]);

features = matrix(iris[:, 2:5]);        # Where is the matrix function?
labels = vector(iris[:, "Species"]);

#
stump = build_stump(labels, features);
print_tree(stump)


# Automated Cascading of Splits

tree = build_tree(labels, features);
print_tree(tree)

# Leaf pruning

length(tree)
pruned = prune_tree(tree, 0.9);
length(pruned)

pruned = prune_tree(tree, 0.6);
length(pruned)
print_tree(pruned)


# Random forest

forest = build_forest(labels, features, 2, 10);
predictions = apply_forest(forest, features);
nfoldCV_forest(labels, features, 2, 10, 3);

# To add extra procesor either restart julia with ./julia -p 2 or
# from within this session add processor addprocs(1) and then load DecisionTree.jl (using DecisionTree)
