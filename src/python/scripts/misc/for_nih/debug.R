
fpath <- '/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_nih/s1/gene_cpm.csv'

df <- read.csv(fpath)

dim(df)
colnames(df)
rownames(df)
df$X

marginal_row_sums <- rowSums(df[,-1])

boxplot(marginal_row_sums)

marginal_column_sums <- colSums(df)
boxplot(marginal_column_sums)