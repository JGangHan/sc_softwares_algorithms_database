

## 1. 整体流程
1. 数据准备：  
   1.1 cellranger 输出的 possorted_genome_bam.bam 文件（samplename/outs/possorted_genome_bam.bam）  
   1.2 Seurat 导出目标 barcode 序列（clean cells，filtered_barcodes.tsv）  
   1.3 参考基因组 gtf 文件（cellranger 参考基因组制作步骤中生成的注释文件 genes.gtf 即可）
2. **Seuurat 导出目标 barcode 序列**  
**因为对原始测序数据比对后的 possorted_genome_bam.bam 文件包含所有的 barcode 序列，其中有的是空白 barcode，有的可能是低质量细胞或多细胞的barcode，如果不指定目标细胞 barcode， velocyto 软件默认会检索所有 barcode，从而浪费大量时间**  
```
library(Seurat)
ob.merge.rmdb = readRDS('ob.merge.rmdb.rds')
# seurat 合并样本过程中会对 barcode 重编码，需要将其恢复，格式如下
# ACAGAAATCAAACGTC-1_1  >>>  ACAGAAATCAAACGTC-1
samplename <- c("E50","E55","E60","E63","E66","E69", "E72","E75","E80")
cb <- ob.merge.rmdb@meta.data
for (i in samplename) {
  # 提取指定样本
  cb_subset <- subset(cb, samplename == i)
  # 修改 row.names，去除 "_N"
  cb_subset <- gsub("_.*", "", row.names(cb_subset))
  
  # 写入文件
  output_filename <- paste0("E:/R working directory/Final and Integration/filtered_barcode/", i, "_filtered_barcodes.tsv")
  write.table(cb_subset, file = output_filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
```
3. velocyto 比对剪切和未剪切 mRNA，最终生成 .loom 文件（每个样本对应一个 loom 文件）  
4. scvelo 进行 RNA 速率分析   
5. 一些有用链接  
[velocyto官网教程](https://github.com/velocyto-team/velocyto-notebooks/blob/master/python/DentateGyrus.ipynb)  
[scVelo在python3.10环境下的安装及使用](https://github.com/PrinceWang2018/scvelo_py3.10?tab=readme-ov-file)  
[基于Seurat UAMP和celltype的RNA Velocity速率分析的全套流程](https://www.jianshu.com/p/c33341b65cad)  



## 2. 软件安装
1. velocyto 安装，
```
conda create --name velocyte   # 默认的 python 3.6 版本
conda activate velocyte
conda install velocyte
# 虽然是一个 python 软件，但可以 linux 环境下直接测试
velocyto -h

# 后续 pysam 报错，重新安装 pysam==0.17 
```
2. scvelo 安装  
这个软件安装起来比较麻烦有很多 python 包不兼容，后期需要自己测试  
```
conda create --name scvelo python=3.9  # python3.6 安不上 sscvelo
conda activate scvelo
conda install bioconda::scvelo  # 默认 scvelo v0.2.5
python3.9 # 注意系统默认 python 路径


import scvelo
# 发生TLS报错，scikit-learn 和 scanpy 版本问题，或者是这两个包与其他依赖包的版本问
# 下列版本兼容  
./pip install scikit-learn==1.5.2
./pip install matplotlib==3.7.1
./pip install scanpy==1.9.3
```

## 3. velocyte 步骤
**因为对原始测序数据比对后的 possorted_genome_bam.bam 文件包含所有的 barcode 序列，其中有的是空白 barcode，有的可能是低质量细胞或多细胞的barcode，如果不指定目标细胞 barcode， velocyto 软件默认会检索所有 barcode，从而浪费大量时间**

### 3.1 测试数据子集
```
# 提取 E50 子集（万分之七）
samtools view -s 0.0007 -@ 16 -b -h /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E50_200150/outs/possorted_genome_bam.bam > /data/hanjiangang/hanjiangang/single_Cell/velocity/test/E50_random_subset.bam
# -b 输出为 bam 文件，-h 保留头部文件
# velocyto 检测剪接和未剪接 mRNA，生成 loom 文件
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E50_filtered_barcodes.tsv -e test -o /data/hanjiangang/hanjiangang/single_Cell/velocity/test /data/hanjiangang/hanjiangang/single_Cell/velocity/test/E50_random_subset.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/test/test.txt &
```
### 3.2 正式代码

```
# -v 输出代码运行信息；-b 目标 barcode 文件；-e 输出文件前缀、输出文件中的细胞barcode前缀和样本信息，最好加上
# E50
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E50_filtered_barcodes.tsv -e E50 -o /data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E50_200150/outs/possorted_genome_bam.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/log/velocyto_E50.txt &
# E55
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E55_filtered_barcodes.tsv -e E55 -o /data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E55_200530/outs/possorted_genome_bam.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/log/velocyto_E55.txt &
# E60
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E60_filtered_barcodes.tsv -e E60 -o /data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E60_200048/outs/possorted_genome_bam.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/log/velocyto_E60.txt &
# E63
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E63_filtered_barcodes.tsv -e E63 -o /data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E63_200136/outs/possorted_genome_bam.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/log/velocyto_E63.txt &
# E66
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E66_filtered_barcodes.tsv -e E66 -o /data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E66_200376/outs/possorted_genome_bam.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/log/velocyto_E66.txt &
# E69
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E69_filtered_barcodes.tsv -e E69 -o /data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E69_200462/outs/possorted_genome_bam.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/log/velocyto_E69.txt &
# E72
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E72_filtered_barcodes.tsv -e E72 -o /data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E72_200254/outs/possorted_genome_bam.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/log/velocyto_E72.txt &
# E75
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E75_filtered_barcodes.tsv -e E75 -o /data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E75_200030/outs/possorted_genome_bam.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/log/velocyto_E75.txt &
# E80
nohup velocyto run -v -b /data/hanjiangang/hanjiangang/single_Cell/velocity/filtered_bc/E80_filtered_barcodes.tsv -e E80 -o /data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom /data/hanjiangang/hanjiangang/single_Cell/data/sc_RNA/E80_200474/outs/possorted_genome_bam.bam /data/hanjiangang/hanjiangang/single_Cell/ref_ovis/genes.gtf > /data/hanjiangang/hanjiangang/single_Cell/velocity/log/velocyto_E80.txt &
```

## 4. scvelo 步骤
**需要将 seurat (R) 对目标细胞的注释信息整合到 anndata (python) 对象中**  
### 4.1 barcode 序列修正，seurat 注释结果， umap 位置结果导出  
**velocyto 分析中会对细胞 barcode 序列进行二次编辑，所以将 seurate barcode 序列对齐到 loom barcode**  
**ACAGAAATCAAACGTC-1_1 (seurat对象) >>>  E50:ACAGAAATCAAACGTC-1 (anndata 对象)**  
```
# R
cb <- ob.merge.rmdb@meta.data
cell_embeddings<-Embeddings(ob.merge.rmdb, reduction = "umap")
colnames(cell_embeddings)

# 1. 提取 umap 信息
cb$UMAP_1 = cell_embeddings[,1]
cb$UMAP_2 = cell_embeddings[,2]
identical(row.names(cb), row.names(cell_embeddings))

# 2. 提取 samplename celltype 等信息
cell_id = cb[, c('samplename', 'seurat_clusters', 'celltype_final', 'UMAP_1', 'UMAP_2')]

# 3. 编辑细胞id
#  id 尾部去掉 _N 并添加 x
cell_id$cell_id<-sapply(row.names(cell_id), 
                        function(x)paste(unlist(strsplit(x, "-"))[1],"x",sep = ""))
# id 头部添加 samplename:
cell_id$cell_id<-paste(cell_id$samplename,
                       cell_id$cell_id,sep = ":")
write.csv(cell_id, file = "cell_id.csv", row.names = FALSE)

# python，读取
cell_id = pd.read_csv("cell_id.csv")
```

### 4.2 数据读取
```
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

# E50  E55  E60  E63  E66  E69  E72  E75  E80 
# 1889 4120 2935 3510 4972 6185 6555 6030 7361 
# 读取数据
sample_E50 = anndata.read_loom("/data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom/E50.loom")
sample_E55 = anndata.read_loom("/data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom/E55.loom")
sample_E60 = anndata.read_loom("/data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom/E60.loom")
sample_E63 = anndata.read_loom("/data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom/E63.loom")
sample_E66 = anndata.read_loom("/data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom/E66.loom")
sample_E69 = anndata.read_loom("/data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom/E69.loom")
sample_E72 = anndata.read_loom("/data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom/E72.loom")
sample_E75 = anndata.read_loom("/data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom/E75.loom")
sample_E80 = anndata.read_loom("/data/hanjiangang/hanjiangang/single_Cell/velocity/velocyto_loom/E80.loom")
```

### 4.3 多样本合并
**因为少数基因（列名）不是唯一出现，同时出现多次，所以合并数据之前需要处理这些基因id**
多样本需要先进行代码测试
```
a= sample_E50
b=sample_E55

# 检查是否存在重复基因
a.var_names.is_unique
# 哪些重复基因
gene_counts = a.var_names.value_counts()
# 筛选出出现次数超过 1 的基因名，表示这些基因名是重复的
duplicate_genes = gene_counts[gene_counts > 1]
duplicate_genes.shape()
# 打印重复基因及其出现次数
print(duplicate_genes)

# 对重复基因重新编号
# a
unique_genes = {}  # 创建一个空字典，用于存储基因名及其出现次数
new_var_names = []  # 创建一个空列表，用于存储新的基因名

for gene in a.var_names:
    # 检查当前基因名是否已经在 unique_genes 字典中
    if gene not in unique_genes:
        unique_genes[gene] = 1  # 如果不在，初始化出现次数为 1
        new_var_names.append(gene)  # 将基因名添加到新列表中
    else:
        unique_genes[gene] += 1  # 如果在，出现次数加 1
        # 将基因名重命名为 "基因名_出现次数" 的格式，并添加到新列表中
        new_var_names.append(f"{gene}_{unique_genes[gene]}")

# 更新 AnnData 对象中的 var_names，将新生成的基因名列表赋值给 a.var_names
a.var_names = new_var_names


# b
unique_genes = {}  # 创建一个空字典，用于存储基因名及其出现次数
new_var_names = []  # 创建一个空列表，用于存储新的基因名
for gene in b.var_names:
    # 检查当前基因名是否已经在 unique_genes 字典中
    if gene not in unique_genes:
        unique_genes[gene] = 1  # 如果不在，初始化出现次数为 1
        new_var_names.append(gene)  # 将基因名添加到新列表中
    else:
        unique_genes[gene] += 1  # 如果在，出现次数加 1
        # 将基因名重命名为 "基因名_出现次数" 的格式，并添加到新列表中
        new_var_names.append(f"{gene}_{unique_genes[gene]}")

# 更新 AnnData 对象中的 var_names，将新生成的基因名列表赋值给 a.var_names
b.var_names = new_var_names

# 检查是否相同
set(a.var_names) == set(b.var_names)
a.var_names.equals(b.var_names)
 
# 合并
samples = [a, b]
c = anndata.concat(samples, join='outer', label='batch', keys=['a', 'b'])
c = anndata.concat(samples, join='inner', label='batch', keys=['a', 'b'])







```




# 9 个对象构建列表
samples = [sample_E50, sample_E55, sample_E60, sample_E63, sample_E66, sample_E69, sample_E72, sample_E75, sample_E80]

# 定义一个函数来处理每个 AnnData 对象的重复基因
def rename_duplicate_genes(adata):
    unique_genes = {}  # 创建一个字典用于存储基因及其出现次数
    new_var_names = []  # 创建一个列表用于存储新的基因名
    # 遍历 var_names，处理重复基因名
    for gene in adata.var_names:
        if gene not in unique_genes:
            unique_genes[gene] = 1
            new_var_names.append(gene)
        else:
            unique_genes[gene] += 1
            new_var_names.append(f"{gene}_{unique_genes[gene]}")
    # 更新 var_names
    adata.var_names = new_var_names

# 循环处理每个 AnnData 对象
for adata in samples:
    rename_duplicate_genes(adata)  # 调用函数处理 var_names

ob_merge = anndata.concat(samples, join='outer', label='batch', keys=['E50', 'E55', 'E60', 'E63', 'E66', 'E69', 'E72', 'E75', 'E80'])
ob_merge.obs.index


# 目标细胞
target = cell_id["cell_id"]
# 筛选合并后的 AnnData 对象，基于 cell_ids_column 中的值
ob_merge = ob_merge[np.isin(ob_merge.obs.index, target)]

ob_merge_index = pd.DataFrame(ob_merge.obs.index) # CellID
ob_merge_index = ob_merge_index.rename(columns={'CellID': 'cell_id'})


# 提取 metadata 中的 'UMAP_1', 'UMAP_2', 和 'cell_id' 列
umap = cell_id[['UMAP_1', 'UMAP_2', 'cell_id']]

# 根据 ob_merge_index 中的 cell_id 顺序对 metadata 进行过滤和排序
# 首先将 metadata 中的 cell_id 过滤，确保它们在 ob_merge_index 中
umap = umap[umap['cell_id'].isin(ob_merge_index['cell_id'])]

# 根据 ob_merge_index 中的 cell_id 进行排序（使用 merge 操作保持顺序）
umap_ordered = ob_merge_index.merge(umap, on='cell_id')
set(umap_ordered[['cell_id']]) == set(ob_merge_index[['cell_id']])  # 检查顺序

# 此时 metadata_ordered 包含按 ob_merge_index 顺序排列的 'UMAP_1' 和 'UMAP_2'
umap_ordered = umap_ordered[['UMAP_1', 'UMAP_2']]

ob_merge.obsm['X_umap'] = umap_ordered.values
ob_merge.obsm



# 添加 celltype 注释信息
celltype = cell_id[['celltype_final', 'cell_id']]
celltype = celltype[celltype['cell_id'].isin(ob_merge_index['cell_id'])]
celltype_ordered = ob_merge_index.merge(celltype, on='cell_id')
set(celltype_ordered[['cell_id']]) == set(ob_merge_index[['cell_id']])  # 检查顺序

# 此时 metadata_ordered 包含按 ob_merge_index 顺序排列的 'UMAP_1' 和 'UMAP_2'
celltype_ordered = celltype_ordered[['celltype_final']]
ob_merge.obs['celltype'] = celltype_ordered.iloc[:, 0].values


# 提取子集
ob_merge.obs['celltype'].value_counts()
example = ob_merge[ob_merge.obs['celltype'].isin(['PSC_C3', 'VSMC_C1', 'VSMC_C2']), :]

## scvelo RNA速率分析
scv.pp.filter_and_normalize(example)
scv.pp.moments(example)
scv.tl.velocity(example, mode = "stochastic")
scv.tl.velocity_graph(example)



```

## scvelo RNA速率分析
scv.pp.filter_and_normalize(ob_merge)
scv.pp.moments(ob_merge)
scv.tl.velocity(ob_merge, mode = "stochastic")
scv.tl.velocity_graph(ob_merge)

# spliced/unspliced的比例
scv.pl.proportions(adata)

# 可视化
scv.pl.velocity_embedding(adata, basis = 'umap')


scv.pl.velocity_embedding_stream(adata, basis = 'umap')



















