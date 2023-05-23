#### RNAseq下游分析流程包【RCB Lab】
##### 1. 得到差异分析基因（DEG）
&ensp;&ensp;**数据格式要求**
1. RNA表达数据，第一列为基因名，建议列名最好设置为gene_name，后续列名为样本编号
2. 分组文件，第一列名为样本编号，列名设置为rowname, 第二列为分组，列名设置为type
3. 注意：即RNA的表达数据的第一行为: gene_name, sample1, sample2, ..., 分组文件的第一列(样本编号)必须要包含于表达数据中的sample1-n。

```R
library(lzRNA)
eset.fname = "../Rpackage_testdata/rna_eset.csv"
group.fname = "../Rpackage_testdata/rna_group.csv"
# eset 处理
eset <- data.table::fread(eset.fname, data.table = F)
if (sum(is.na(eset)) > 0) {
  eset <- na.omit(eset)
  warning("Your exprs data have NA value.\n...Delete NA")
}
eset <- DEG_preEset.dup(df = eset, index = "gene_name")
# group 处理
row.select = NULL
group <- readr::read_csv(group.fname, show_col_types = F)
if( !is.null(row.select) ) {
  group <- group[unlist(strsplit(row.select, split = "")),]
}
group <- as.data.frame(group, check.names = F)
group <- column_to_rownames(group, var = 'rowname')
group$type <- factor(group$type)

# 2. 差异分析
cat("分析中...")
eset.group.list <- list(eset = eset <- eset[, rownames(group)] ,
                        group = group, f_mark = '1')
diffan.d <- DEG_edgeR(eset.group.list)
```