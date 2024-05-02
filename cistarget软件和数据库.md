## [cistarget数据库](https://resources.aertslab.org/cistarget/)
## [Rcistarget](https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html)
## [pycistarget](https://pycistarget.readthedocs.io/en/latest/tools.html#)



## cisTarget 数据库
主要是一个利用已知转录因子motif序列信息，针对每一个motif（数量超过20k），利用**cluster-buster软件**检索人/小鼠/果蝇所有基因TSS位点上游500bp/上下游5kb/上下游10kb，按照所有基因检索序列与该motif序列结合的可能性的大小（cbust软件输出的CRM scores）进行排序，生成 genes ranking-motif matrix文件，用于基因的motif富集分析。还有一种情况是在检索基因TSS上下游序列过程中，可能会有序列片段与motif基序匹配，因此生成 regions ranking-motif matrix文件，可以用于序列的motif富集分析
![image](https://github.com/JGangHan/Software-list/assets/75400599/777fbf17-31f9-4b50-812d-685fe20582ce)

## 构建专属物种的cisTarget 数据库




## Rcistarget 分析流程
1. 输入文件
   * 目标基因集 gene set
   * gene-motif ranking 文件
   * motif 注释文件：用于注释motif序列对应的转录因子
2. 什么是 gene-motif ranking 文件
   * 是一个利用已知转录因子motif序列信息，针对每一个motif（数量超过20k），检索所有基因TSS位点上游500bp/上下游5kb/上下游10kb，按照所有基因检索序列与该motif序列结合的可能性的大小进行排序
   * 目前只有人、小鼠和果蝇的
   * 可以自己构建专属物种的gene-motif ranking 文件
3. 计算AUC score：在每一个motif中，选择top 5%可能结合的基因，构建目标基因集-top genes的恢复曲线，从而计算目标基因集在每个motif中的AUC score
4. 所以输出矩阵为：行表示基因集，列表示基序，矩阵中的每个值代表对应基因集和基序的AUC值。下图为单个基因集AUC score的分布特征，为长尾分布
  ![003b4647d93490d1095c178451478f79_wOl1Ts64OJFcwAAAABJRU5ErkJggg==](https://github.com/JGangHan/Software-list/assets/75400599/98138926-4805-4d23-a6b6-eb037abfb909)
5. 对AUC score进行矫正得到Normalized Enrichment Score (NES)，设定阈值（默认为3）筛选最显著富集到的motif
6. 对motif进行注释，得到转录因子信息
7. 需要注意**这里是对一个gene list中的多个基因进行motif富集分析，并不意味着所有的基因都能结合富集到的转录因子**，因此可以具体分析是某个或某些基因能结合到富集到的转录因子
下图为：显著富集到的motif对基因集中的基因恢复情况
![image](https://github.com/JGangHan/Software-list/assets/75400599/261217f3-e6be-47f6-9459-51c18c0c0589)

### gene-motif ranking 文件
**让人困惑的地方**：
1. 首先就是名字，文件的本身内容是 motif*gene ranking矩阵，这里却叫gene-motif ranking 矩阵
2. 在针对单细胞分析的AUCell流程中，为gene ranking * cell矩阵，这个文件内容却是 motif * gene ranking矩阵，两者行列完全颠倒，在同时理解这两个矩阵的时候，让人困惑
3. **这两个矩阵的实质上都是对基因的排序，只不过AUCell中是根据所有基因在单个细胞中的基因表达水平，对所有基因进行排序；而Rcistarget中是根据所有基因与单个motif结合可能性的大小，对基因进行排序。之后构建目标基因集-top genes（5%）恢复曲线占比情况表示富集程度**
![image](https://github.com/JGangHan/Software-list/assets/75400599/3a4d0477-3bef-4a5c-9187-3acb4d54bc5a)

### motif annotation 文件
![image](https://github.com/JGangHan/Software-list/assets/75400599/669e2a3f-46a2-410c-8a08-2dc08b984ec1)



## pycistarget
与Rcistarget不同，Rcistarget更倾向于对基因集进行motif富集分析，而pycistarget更倾向于对筛选到的基因序列/区域进行motif富集分析
















