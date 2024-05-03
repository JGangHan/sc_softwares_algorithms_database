## [cistarget数据库](https://resources.aertslab.org/cistarget/)
## [Rcistarget](https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html)
## [pycistarget](https://pycistarget.readthedocs.io/en/latest/tools.html#)

## pycistarget tips
* **与Rcistarget不同，Rcistarget 更倾向于对基因集进行motif富集分析，而pycistarget更倾向于对筛选到的基因序列/区域进行motif富集分析**
* 在Rcistarget中可以对不同的基因集进行motif富集分析，在Pycistarget中可以
  * 对多个 regions lists 进行 motif 富集分析，
  * 并可进一步比较不同 regions lists 的富集结果，即 Differentially Enriched Motifs (DEM) 分析
* **没有弄清楚以下问题**
  1. cisTarget 数据库仅提供了人/小鼠/果蝇这三种模式生物的 regions vs motifs ranking matrix，当对其他物种进行 motif 富集分析时，确实需要构建自定义 cistarget 数据库，在后续分析过程中由于 regions 是完全对应关系，所以不需要担心自定义 cistarget 数据库和目标 regions 不一致情况
  2. 但对人/小鼠/果蝇这三种模式生物进行分析时，目标 regions 绝对与[用于构建 cistarget 数据库所使用的 regions](https://resources.aertslab.org/cistarget/regions/) 不一样，那是否仍然需构建自定义数据库呢？
     * 如果不许需要的话，是不是证明在后续对目标 regions 进行 motif 富集分析时，pycistarget 软件会对 目标 regions 和 数据库 regions 进行比较和处理？

       
## pycistarget motif 富集分析流程
1. 准备 cisTarget 数据库，即 regions vs motifs ranking matrix
2. 准备 regions sets：可以包含多个 regions list
3. 设定阈值，进行 motif 富集分析，其大概原理应该是与 Rcistarget 类似，**观察目标 regions list 中的 regions 在 regions vs motifs ranking matrix 中各 motif 基于对应的 top regions 中的恢复程度，最后每个 motif 基于都会有一个对应该 regions list 的 AUC score**
4. 输出符合阈值的 motif 基序还有对应的 motif to tfs 注释信息
5. 可以观察指定的 motif 基序（必须是富集到的）对应的潜在结合 regions


## pycistarget Differentially Enriched Motifs (DEM) 分析流程
**比较两个 regions list 所富集到的 motifs 差异**
1. 






