julia summary/collect.jl
Rscript summary/plots-tables.R
../data/sims

# find ../data/sims/ -name '*.png' | cpio -pdm  ../../CautiousBootstrap/papers/main-paper/figures/cautious-learning-fixed/
cp -r ../plots/sims ../../CautiousBootstrap/papers/main-paper/figures/
cp -r ../plots/ICUadmissions ../../CautiousBootstrap/papers/main-paper/figures/
