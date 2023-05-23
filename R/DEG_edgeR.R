#' RNAseq differential analysis
#'
#' @description The function is designed to perform differential analysis by edgeR.
#'
#' @details Following data is requried: express, group.
#'
#' @param exprset.group x is a 3-object list: list(eset=exprset, group=phen, f_mark="")
#' @param coef default is 2, meanig the 2nd fator vs 1st factor of group.
#' @param pval default is 0.05
#' @param fdr default is 0.1
#' @param logfc default is 1
#'
#' @return a list: list(qlf = qlf, resdf = et.norm, deg = etSig)
#' qlf: mid data of edgeR, et.norm: all result of DEG, sig result of DEG.
#'
#' @export
#'
#' @importFrom stats median model.matrix
#' @import edgeR
#' @import dplyr
#'
#' @examples # Run: diffan <- DEG_edgeR(exprset.group) # pval=0.05, fdr=0.1, logfc=1
#' # exprset.group  由下列三个数据打包而成，示例：list(eset=exprset, group=phen, f_mark="")
#' # *1 exprset.group$eset： 矩阵数据格式(数值型，整型)
#' #       | row1 | row2 | row3  | row4
#' # gene1 |  34  |  23  |  56   |  23
#' # gene2 |  35  |  23  |  12   |  23
#' # gene3 |  12  |  78  |  78   |  78
#' # *2 exprset.group$group： 分组数据格式：建议只有一列数据
#' #    # 需要组的行名=表达谱的列名 rownames(group) == colname(eset)
#' #            type
#' # rowname1 | tumor
#' # rowname2 | tumor
#' # rowname3 | normal
#' # rowname4 | normal
#' # *3 exprset.group$f_mark：list标记项，可以为空，但是如果是批量差异分析，建议必须要有此参数，方便后续对个结果进行保存和导出
DEG_edgeR <- function(exprset.group, coef = 2,
                      pval=0.05, fdr=0.1, logfc=1) {
  # 此函数为edegR差异分析，返回两个表格组成的list：
  #   1.全部的分析表格 + 2.差异表格(P<0.05&FDR<0.1&abs(et$log2FC)>1)
  # exprset.group要求为list (主要为了方便lapply函数做批量差异分析)

  # 1.0 count，group
  cat('\n', exprset.group$f_mark, ' ================\n')
  exprset <- exprset.group$eset
  pheno <- exprset.group$group
  exprset <- exprset[, rownames(pheno)]

  # 1.1 构建DEGList
  y <- DGEList(counts = exprset, group = pheno[,1])
  keep <- filterByExpr(y)
  #keep <- filterByExpr(y, min.count = 10, min.total.count = 15,
  #                     large.n = 10, min.prop = 0.7)
  y <- y[keep, , keep.lib.sizes=FALSE]
  # 进行TMM标准化
  y.norm <- calcNormFactors(y, method = "TMM")
  # 获取cpm TMM值
  cpm.tmm <- cpm(y.norm, log = F)
  cpm.tmm <- log2(cpm.tmm + 1)

  # 1.2 计算离散因子
  y <- calcNormFactors(y)
  # design.form <- as.formula(paste("~", "group"))
  # # y$samples$group <- relevel(y$samples$group, ref="C")
  design <- model.matrix(~group, data =  y$samples)
  y <- estimateDisp(y, design)
  #To perform quasi-likelihood F-tests:
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = coef)

  # 获取结果
  # 获取排名靠前的基因，这里设置n=80000是为了输出所以基因
  et <- topTags(qlf, n = 80000)
  et <- as.data.frame(et) # 转换为数据框类型
  et.norm <- merge(et, cpm.tmm, by = 0)
  et.norm <- column_to_rownames(et.norm, var = "Row.names")
  et.norm$Gene <- rownames(et.norm)
  et.norm <- et.norm %>% dplyr::rename(log2FC=logFC) %>%
    dplyr::arrange(desc(log2FC), PValue) %>%
    dplyr::select(Gene, log2FC, PValue, FDR, everything())
  # 差异基因筛选
  etSig <- et.norm[which(et.norm$PValue < pval &
                           et.norm$FDR < fdr & abs(et.norm$log2FC) > logfc),]
  return(list(qlf = qlf, resdf = et.norm, deg = etSig))
}



#' removes duplicate row from a dataframe
#'
#' @description removes duplicate row from a dataframe, retains the row with the higher mean.
#'
#' @param df dataframe the first of the dataframe is the rowindex and recommend to set name as gene_name.
#' @param index character default is "gene_name", if the fist column of df is not gene_name, please set the column name.
#'
#' @return dataframe a express eset which rowname is the unduplicated index
#' @export
#'
#' @import dplyr
#' @import tibble
#'
#' @examples # eset <- DEG_preEset.dup(df = eset) # index = "gene_name"
#' # df:
#' # |gene_name | row1 | row2 | row3  | row4
#' # | gene1    |  34  |  23  |  56   |  23
#' # | gene2    |  35  |  23  |  12   |  23
#' # | gene3    |  12  |  78  |  78   |  78
DEG_preEset.dup <- function(df, index="gene_name") {
  if (sum(duplicated(df[, 1])) > 0) {
    df$median <- apply(df[, 2:ncol(df)], 1, median)
    df$mean <- rowMeans(df[,2:ncol(df)])
    df <- df %>% arrange(desc(median), desc(mean))
    df <- df[!duplicated(df[, 1]), ]
    rownames(df) <- NULL
    df <- column_to_rownames(df, var = index)
    df$median <- NULL;df$mean <- NULL
    return(df)
  } else {
    rownames(df) <- NULL
    df <- column_to_rownames(df, var = index)
    return(df)
  }
}
