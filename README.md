# 软件/包
## SCENIC R/Python 基因表达网络分析
基于scRNA-seq基因表达数据和人/小鼠转录因子列表，首先针对转录因子进行基因共表达分析，之后利用TFs motif序列和基因上游调控序列对TFs共表达基因进行过滤，最后在每个细胞中计算AUC指数。这个软件本质是**基因共表达分析**，在其基础之上加上了motif序列过滤和AUC富集指数计算

## SCENIC+ Python

## GENIE3 
基因共表达分析


## [Cluster-Buster (cbust)](https://bu.wenglab.org/cluster-buster/index.html)，在目标DNA片段中检索可能存在的motif基序或motif cluster
利用多个motif的位置权重矩阵，cbust可以在多个DNA片段序列中检索可能存在的motif cluster及其对应的motifs，并输出motif/motifs cluster 区间，DNA链，评分等多项评估信息


## [RcisTarget](https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html); [pycistarget](https://pycistarget.readthedocs.io/en/latest/tools.html#); [cistarget数据库](https://resources.aertslab.org/cistarget/)
* **cisTarget** 只是一个数据库，对 Rcistarget 或 Pcistarget motif 富集分析提供背景数据库
  利用已知转录因子motif序列信息，针对每一个motif（数量超过20k），检索人/小鼠/果蝇所有基因TSS位点上游500bp/上下游5kb/上下游10kb，按照所有基因检索序列与该motif序列结合的可能性的大小进行排序，生成 genes ranking-motif matrix文件
* **Rcistarget**将单个或多个基因集，与 genes ranking-motif matrix 中的每一个motif所包含的所有经过排序的基因构建恢复曲线（revovery curve），从而计算基因集在每一个motif的AUC score，之后矫正为normalized enrichemnt score （NES），将筛选到的显著富集motif基序注释为转录因子
* 






 i-cisTarget and iRegulon

## [AUCell，Area Under the Curve（AUC）](https://www.bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html) 
评估基因集在单个细胞中的活性/富集程度（也就是AUC score），例如某些个调控子regulons的活性，从而生成一个AUC socres*cells数据集
注：调控子regulon：一个regulon指的某一个转录因子及其假定靶基因的集合



## [CoveringPackage, 筛选可用于注释细胞的 gene list](https://github.com/lanlanji/CoveringPackage/tree/master)
将 raw count matrix 粗暴的转为二进制矩阵，基于已有的细胞注释，通过 CellCover 算法筛选对每个细胞类型最优的 markers panel（也就是gene list），然后验证该 markers panel 在每个细胞类群中的覆盖度 covering rate，理论上是markers panel 在目标细胞类型中的最高，从而实现用 gene list 来注释细胞类型

## [MAST, 单细胞基因差异分析](https://lishensuo.github.io/posts/bioinfo/030%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--mast%E5%B7%AE%E5%BC%82%E5%9F%BA%E5%9B%A0%E5%88%86%E6%9E%90/)
综合考虑多种协变量因素，例如离群细胞，基因在细胞中的表达比例等



# 算法/函数
## Jensen-Shannon divergency and distance, JS散度（JSD）和JS距离，可以用于计算两个向量空间距离或相似性，或基因表达/活性特异性    
1. Jensen-Shannon 散度是一种流行的度量，用于衡量两个概率分布之间的相似性，其特征是对称性和有界性。
   *对称性指的是JSD(P||Q) = JSD(Q||P)，也就是无论如何交换两个分布的位置，计算结果相同；
   *有界性指的是JSD值位于0-1之间，JSD=0，表示两个分布完全相同，JSD=1表示两个分布完全不同
2. Jensen-Shannon距离Jensen-Shannon 散度的平方根，使其成为一个真正的距离度量，其值也在 0 到 1 之间
3. 在python中可以通过from scipy.spatial.distance import jensenshannon调用jensenshannon()函数
4. Scenic软件中在计算[Regulon Specificty Scores (RSS)指数](https://github.com/aertslab/pySCENIC/blob/master/src/pyscenic/rss.py)调用了该参数，其具体实现思路是
   * 使用AUC matrix（可以理解为多个向量的集合）
   * 将细胞标签矩阵转换为二进制矩阵（0/1矩阵），在分析过程中会提取每一列的细胞注释结果，属于目标细胞类型为1，不属于为0 
   * 用 “1.0 - jensenshannon(aucs / aucs.sum(), labels / labels.sum())”表示RSS，也就是说RSS值越接近1，表明在目标细胞类型中特异性越强


## colRanks()函数，计算矩阵中每一列元素的排名
AUCell软件通过 colRanks() 函数先将表达矩阵（genes * cells）转换为排序矩阵(ranking * cells)，进而计算 AUC score*cells矩阵


## 阈值计算策略，AUCell/R/priv_auc.assignmnetThreshold_v6.R
AUCell包的源码中提供了自动计算阈值的方法，基于多种统计方法来计算和选择最优阈值，如果以后有涉及阈值选择的问题可以参考该方法


## Adjusted Rand Index, ARI, 衡量两个矩阵聚类结果相似性
ARI 是一种用于衡量**两个数据聚类结果相似性**的统计度量，它是 Rand Index 的一个改进版本。Rand Index 本身用于评估两种数据聚类之间的相似度，但它没有考虑到随机分类的影响。调整后的 Rand Index 通过对原始 Rand Index 进行调整，考虑到了随机分类可能产生的影响，使得结果更加准确和可靠。**ARI的取值范围从-1到1，可以通过 R 的 mclust/cluster 包，或 python 的 scikit-learn 包计算**。
* 当ARI为1时，表示两个聚类完全一致。
* 当ARI为0时，表示聚类结果与随机聚类没有区别，即聚类效果与随机效果相同。
* 当ARI为负值时，表示聚类结果比随机聚类更差。
例如：CellCover 中比较了原始矩阵和由原始矩阵转化的二进制矩阵，使用不同基因数量的聚类结果

![image](https://github.com/JGangHan/software_information/assets/75400599/d16d7cb1-d830-4afc-9dae-298442247086)



# 数据库
## [CZ Cell*Gene discovery](https://cellxgene.cziscience.com/docs/01__CellxGene) 人/小鼠单细胞数据库
![image](https://github.com/JGangHan/sc_softwares_algorithms_database/assets/75400599/c0818727-321f-4cea-a214-86cabce9e3e9)

1. 主要是作为一个各种单细胞测序数据的数据集合，非常庞大，可免费下载 .h5ad 格式和 .rds 格式处理后的单细胞数据
2. 样本信息包含组织、年龄、生理状态、细胞类型等
3. 还提供各种网页工具探究基因表达，数据下载工具，注释工具
