## [cistarget数据库](https://resources.aertslab.org/cistarget/)
## [Rcistarget](https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html)
## [pycistarget](https://pycistarget.readthedocs.io/en/latest/tools.html#)



## cisTarget 数据库
* 主要是一个利用已知转录因子motif序列信息，针对每一个motif（数量超过20k），利用**cluster-buster软件**检索人/小鼠/果蝇所有基因TSS位点上游500bp/上下游5kb/上下游10kb 区间范围
* 按照所有基因检索序列与该motif序列结合的可能性的大小（**cbust软件输出的CRM scores**）进行排序，
* 生成 **regions ranking-motif matrix文件**，该文件可以//通过Pycistarget软件用于基因序列motif富集分析**；
* 将属于同一基因的regions合并，即可生成 **genes ranking-motif matrix 文件**，该文件可通过**Rcistarget包用于基因的motif富集分析**。
* **其核心是region vs motif score matrix，和基于得分进一步排序得到的 region vs motif raning matrix**
  
**以人为例：**
* 检索生成的score文件：hg38_screen_v10_clust.regions_vs_motifs.scores.feather
* 根据score文件对区间进行排序：hg38_screen_v10_clust.regions_vs_motifs.rankings.feather

**数据库目录结构：**
![image](https://github.com/JGangHan/Software-list/assets/75400599/777fbf17-31f9-4b50-812d-685fe20582ce)

**genes vs motif.ranking文件格式**
![image](https://github.com/JGangHan/software_information/assets/75400599/8422f63b-9c26-4f97-815a-b9d2b364d182)



## [构建专属物种的cisTarget数据库](https://github.com/aertslab/create_cisTarget_databases?tab=readme-ov-file#create_cistarget_motif_databasespy)
### 1. 主要实现下方六种功能
**主要不太理解这个track是什么意思，可能是涉及‘bigWig files of TF ChIP-seq data’的一种文件格式，可能是用track（类似于motif的一种特殊序列，或潜在的新的motif基序），由于之前没接触过，这里跳过**
* **create_cistarget_motif_databases.py：针对motif构建cistarge database**
* create_cistarget_track_databases.py：
* combine_partial_motifs_or_tracks_vs_regions_or_genes_scores_cistarget_dbs.py：由于内存限制或时间限制拆分构建
* combine_partial_regions_or_genes_vs_motifs_or_tracks_scores_cistarget_dbs.py：同上
* convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py：将 score matrix 转成 ranking matrix
* create_cross_species_motifs_rankings_db.py：构建跨物种cistarge database

### 2.在构建绵羊cisTarget数据库过程中仅使用了‘create_cistarget_motif_databases.py’
**Tips**
1. 因为是需要用SCENIC做scRNA-seq+scATAC-seq联合分析，所以是通过peak区域构建的regions vs motif ranking matrix
2. 可以使用其他类型的目标DNA区域
3. 如果是构建genes vs motif ranking matrix，类似cisTarget示例数据集，应该是需要预先提取基因TSS上下游目标区域DNA序列（因为没做过，所以不太确定）

**大概流程**
1. 如果 scATAC-seq 数据分析过程用 pycisTopic 分析，会直接对 peaks 生成 consensus_regions.bed 文件
2. 根据 consensus_regions.bed 文件提取序列 consensus_regions.fa
3. motif 基序数据集 https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
4. 染色体格式 chr1
5. 确认好输入文件格式、文件路径和软件路径之后就可以直接运行
```
/home/hanjiangang/anaconda3/envs/create_cistarget_databases/bin/python3.10 \
/data/hanjiangang/sc/pycistarget/create_cisTarget_databases/create_cistarget_motif_databases.py \
-f $consensdir \
-M $cbdir \
-m $motif_list \
-o $outdir/$tag \
-t 50 \
-c $cbpath \
-l
```
5. 输出三个文件
   * cluster_V10_DPCL_sheep.motifs_vs_regions.scores.feather
   * **cluster_V10_DPCL_sheep.regions_vs_motifs.rankings.feather** **该文件直接用于对regions进行motif富集分析**
   * cluster_V10_DPCL_sheep.regions_vs_motifs.scores.feather













