## [cistarget数据库](https://resources.aertslab.org/cistarget/)
## [Rcistarget](https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html)
## [pycistarget](https://pycistarget.readthedocs.io/en/latest/tools.html#)



## pycistarget
* **与Rcistarget不同，Rcistarget 更倾向于对基因集进行motif富集分析，而pycistarget更倾向于对筛选到的基因序列/区域进行motif富集分析**
* 在Rcistarget中可以对不同的基因集进行motif富集分析，在Pycistarget中可以对不同的序列集（regions sets）进行 motif 富集分析，并可进一步进行 Differentially Enriched Motifs (DEM) 分析

**分析流程**
1. 准备 cisTarget 数据库，即 regions vs motifs ranking matrix
2. 准备 regions sets：可以包含多个 regions list
3. 设定阈值，进行 motif 富集分析，其大概原理应该是与 Rcistarget 类似，**观察目标 regions list 中的 regions 在 regions vs motifs ranking matrix 中各 motif 基于对应的 top regions 中的恢复程度，最后每个 motif 基于都会有一个对应该 regions list 的 AUC score**
4. 输出符合阈值的 motif 基序还有对应的 motif to tfs 注释信息
5. 可以观察指定的 motif 基序（必须是富集到的）对应的潜在结合 regions










