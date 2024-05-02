## cistarget数据库（https://resources.aertslab.org/cistarget/）
## [Rcistarget](https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html)
## [pycistarget](https://pycistarget.readthedocs.io/en/latest/tools.html#)



## Rcistarget 分析流程
1. 输入文件
   * 目标基因集 gene set
   * gene-motif ranking 文件
   * motif 注释文件：用于注释motif序列对应的转录因子
2. 什么是 gene-motif ranking 文件
   * 是一个利用已知转录因子motif序列信息，针对每一个motif（数量超过20k），检索所有基因TSS位点上游500bp/上下游5kb/上下游10kb，按照所有基因检索序列与该motif序列结合的可能性的大小进行排序
   * 目前只有人、小鼠和果蝇的
   * 可以自己构建专属物种的gene-motif ranking 文件
3. 在每一个motif中，




## gene-motif ranking 文件
**让人困惑的地方**：
1. 首先就是名字，文件的本身内容是 motif*gene ranking矩阵，这里却叫gene-motif ranking 矩阵
2. 在针对单细胞分析的AUCell流程中，为gene ranking * cell矩阵，这个文件内容却是 motif * gene ranking矩阵，两者行列完全颠倒，在同时理解这两个矩阵的时候，让人困惑
3. **这两个矩阵的实质上都是对基因的排序，只不过AUCell中是根据所有基因在单个细胞中的基因表达水平，对所有基因进行排序；而Rcistarget中是根据所有基因与单个motif结合可能性的大小，对基因进行排序。之后构建目标基因集-top genes（5%）恢复曲线占比情况表示富集程度**
![image](https://github.com/JGangHan/Software-list/assets/75400599/3a4d0477-3bef-4a5c-9187-3acb4d54bc5a)

## motif annotation 文件
![image](https://github.com/JGangHan/Software-list/assets/75400599/669e2a3f-46a2-410c-8a08-2dc08b984ec1)









