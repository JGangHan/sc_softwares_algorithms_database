

## 1. 整体流程
1. 数据准备：  
   1.1 cellranger 输出的 possorted_genome_bam.bam 文件（samplename/outs/possorted_genome_bam.bam）  
   1.2 Seurat 导出目标 barcode 序列（clean cells，filtered_barcodes.tsv）  
   1.3 参考基因组 gtf 文件（cellranger 参考基因组制作步骤中生成的注释文件 genes.gtf 即可）
2. **Seuurat 导出目标 barcode 序列**  
**因为对原始测序数据比对后的 possorted_genome_bam.bam 文件包含所有的 barcode 序列，其中有超级多空白 barcode，有的可能是低质量细胞或多细胞的barcode，如果不指定目标细胞 barcode， velocyto 软件默认会检索所有 barcode，从而浪费大量时间**  
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
[scVelo 官网](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html#)


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

## 4. 数据预处理
**为了保证最终结果的一致性，比如提取指定细胞类型、保证各细胞空间分布，必须将 seurat (R) 对目标细胞的注释信息整合到 anndata (python) 对象中**  
### 4.1 barcode 序列修正，seurat 注释结果， umap 位置结果导出  
**velocyto 分析中会对细胞 barcode 序列进行二次处理，所以后续需要将 seurate barcode 序列对齐到 loom barcode**  
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

### 4.3 多样本合并，测试数据
**因为少数基因（列名）不是唯一出现，同时出现多次，所以合并数据之前需要处理这些基因id**
多样本最好先进行代码测试，多样本合并总是报错，最后发现存在多个相同的 gene symbol
```
a= sample_E50
b=sample_E55

# 1. 检查是否存在重复基因
a.var_names.is_unique
# 哪些重复基因
gene_counts = a.var_names.value_counts()
# 筛选出出现次数超过 1 的基因名，表示这些基因名是重复的
duplicate_genes = gene_counts[gene_counts > 1]
duplicate_genes.shape()
# 打印重复基因及其出现次数
print(duplicate_genes)

# 2. a: 对重复基因重新编号
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

# 3. b: 对重复基因重新编号
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

# 4. 检查ab基因排序是否相同
set(a.var_names) == set(b.var_names)
a.var_names.equals(b.var_names)
 
# 5. 成功合并
samples = [a, b]
c = anndata.concat(samples, join='outer', label='batch', keys=['a', 'b'])
c = anndata.concat(samples, join='inner', label='batch', keys=['a', 'b'])
```


### 4.4 多样本合并，9个正式样本
```
# 9 个对象构建列表
samples = [sample_E50, sample_E55, sample_E60, sample_E63, sample_E66, sample_E69, sample_E72, sample_E75, sample_E80]

# 定义一个函数来处理每个样本数据（AnnData 对象）的重复基因
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

# 合并
ob_merge = anndata.concat(samples, join='outer', label='batch', keys=['E50', 'E55', 'E60', 'E63', 'E66', 'E69', 'E72', 'E75', 'E80'])
ob_merge.obs.index  # 20477*43557，与 seurat 中的数据维度相同
ob_merge.write('./ob_merge_first.h5ad')
a = anndata.read('./ob_merge_first.h5ad')
```


### 4.5 seurat 细胞信息整合
1. 先检查 seurat 细胞和 anndata 细胞是否一致并对其
2. 添加 umap 结果（如果涉及不同细胞子集的重新构图，为保证异质性，后续需要在子集中更新 umap 数据）
3. 添加细胞注释结果
```
# 1. 检查数据是否一致
# seurat 目标细胞（4.1步骤）
target = cell_id["cell_id"]
ob_merge = ob_merge[np.isin(ob_merge.obs.index, target)]

# 2. 提取 UMAP 信息，并对其细胞排列顺序
ob_merge_index = pd.DataFrame(ob_merge.obs.index) # CellID
ob_merge_index = ob_merge_index.rename(columns={'CellID': 'cell_id'})
# umap 子集
umap = cell_id[['UMAP_1', 'UMAP_2', 'cell_id']]
# 检查是否一致
umap = umap[umap['cell_id'].isin(ob_merge_index['cell_id'])]
# 对齐顺序
umap_ordered = ob_merge_index.merge(umap, on='cell_id')
set(umap_ordered[['cell_id']]) == set(ob_merge_index[['cell_id']])  # 检查顺序
# 仅保留 UMAP_1 和 UMAP_2 两列
umap_ordered = umap_ordered[['UMAP_1', 'UMAP_2']]

# 3. 添加 UMAP 信息
ob_merge.obsm['X_umap'] = umap_ordered.values
ob_merge.obsm

# 4. 添加 celltype 注释信息
celltype = cell_id[['celltype_final', 'cell_id']]
celltype = celltype[celltype['cell_id'].isin(ob_merge_index['cell_id'])]
celltype_ordered = ob_merge_index.merge(celltype, on='cell_id')  # 排序
set(celltype_ordered['cell_id']) == set(ob_merge_index['cell_id']) # 检查顺序
set(ob_merge_index['cell_id']) == set(ob_merge.obs.index) # 排序

# 仅保留一列
celltype_ordered = celltype_ordered[['celltype_final']]

# 添加细胞注释结果
ob_merge.obs['celltype'] = celltype_ordered.iloc[:, 0].values

# 5. 数据保存
ob_merge.write('./ob_merge_second.h5ad')
```


### 5. scvelo RNA速率分析，scv.tl.velocity_graph 代码修正
1. 因为服务器配置问题，在 python 环境下无法使用数据并行处理，例如 multiprocessing, parallel。而在 scvelo RNA 速率分析过程中，尤其是 scv.tl.velocity_graph 步骤，会自动调用数据并行处理相关包，所以必须将**velocity_graph 中并行处理相关代码修改为简单命令行**
2. 若不修改代码直接停止运行，或报错**ConnectionResetError：[Errno 104] Connection reset by peer**
3. 修改后代码如下，保存至原路径
```
# 初始代码地址：https://github.com/theislab/scvelo/blob/bef8f1c2044d65a6c6b11e2f519ffafbd7dc2b34/scvelo/tools/velocity_graph.py#L265
import os
import numpy as np
from scipy.sparse import coo_matrix, issparse
from scvelo import logging as logg
from scvelo import settings
from scvelo.core import l2_norm  # 删除了两个需要导入模块
from scvelo.preprocessing.moments import get_moments
from scvelo.preprocessing.neighbors import (
    get_n_neighs,
    get_neighs,
    neighbors,
    pca,
    verify_neighbors,
)
from .utils import cosine_correlation, get_indices, get_iterative_indices
from .velocity import velocity


def vals_to_csr(vals, rows, cols, shape, split_negative=False):
    graph = coo_matrix((vals, (rows, cols)), shape=shape)

    if split_negative:
        graph_neg = graph.copy()

        graph.data = np.clip(graph.data, 0, 1)
        graph_neg.data = np.clip(graph_neg.data, -1, 0)

        graph.eliminate_zeros()
        graph_neg.eliminate_zeros()

        return graph.tocsr(), graph_neg.tocsr()

    else:
        return graph.tocsr()


class VelocityGraph:
    def __init__(
        self,
        adata,
        vkey="velocity",
        xkey="Ms",
        tkey=None,
        basis=None,
        n_neighbors=None,
        sqrt_transform=None,
        n_recurse_neighbors=None,
        random_neighbors_at_max=None,
        gene_subset=None,
        approx=None,
        report=False,
        compute_uncertainties=None,
        mode_neighbors="distances",
    ):

        subset = np.ones(adata.n_vars, bool)
        if gene_subset is not None:
            var_names_subset = adata.var_names.isin(gene_subset)
            subset &= var_names_subset if len(var_names_subset) > 0 else gene_subset
        elif f"{vkey}_genes" in adata.var.keys():
            subset &= np.array(adata.var[f"{vkey}_genes"].values, dtype=bool)

        xkey = xkey if xkey in adata.layers.keys() else "spliced"

        X = np.array(
            adata.layers[xkey].A[:, subset]
            if issparse(adata.layers[xkey])
            else adata.layers[xkey][:, subset]
        )
        V = np.array(
            adata.layers[vkey].A[:, subset]
            if issparse(adata.layers[vkey])
            else adata.layers[vkey][:, subset]
        )

        nans = np.isnan(np.sum(V, axis=0))
        if np.any(nans):
            X = X[:, ~nans]
            V = V[:, ~nans]

        if approx is True and X.shape[1] > 100:
            X_pca, PCs, _, _ = pca(X, n_comps=30, svd_solver="arpack", return_info=True)
            self.X = np.array(X_pca, dtype=np.float32)
            self.V = (V - V.mean(0)).dot(PCs.T)
            self.V[V.sum(1) == 0] = 0
        else:
            self.X = np.array(X, dtype=np.float32)
            self.V = np.array(V, dtype=np.float32)
        self.V_raw = np.array(self.V)

        self.sqrt_transform = sqrt_transform
        uns_key = f"{vkey}_params"
        if self.sqrt_transform is None:
            if uns_key in adata.uns.keys() and "mode" in adata.uns[uns_key]:
                self.sqrt_transform = adata.uns[uns_key]["mode"] == "stochastic"
        if self.sqrt_transform:
            self.V = np.sqrt(np.abs(self.V)) * np.sign(self.V)
        self.V -= np.nanmean(self.V, axis=1)[:, None]

        self.n_recurse_neighbors = n_recurse_neighbors
        if self.n_recurse_neighbors is None:
            if n_neighbors is not None or mode_neighbors == "connectivities":
                self.n_recurse_neighbors = 1
            else:
                self.n_recurse_neighbors = 2

        if "neighbors" not in adata.uns.keys():
            neighbors(adata)
        if np.min((get_neighs(adata, "distances") > 0).sum(1).A1) == 0:
            raise ValueError(
                "Your neighbor graph seems to be corrupted. "
                "Consider recomputing via pp.neighbors."
            )
        if n_neighbors is None or n_neighbors <= get_n_neighs(adata):
            self.indices = get_indices(
                dist=get_neighs(adata, "distances"),
                n_neighbors=n_neighbors,
                mode_neighbors=mode_neighbors,
            )[0]
        else:
            if basis is None:
                basis_keys = ["X_pca", "X_tsne", "X_umap"]
                basis = [key for key in basis_keys if key in adata.obsm.keys()][-1]
            elif f"X_{basis}" in adata.obsm.keys():
                basis = f"X_{basis}"

            if isinstance(approx, str) and approx in adata.obsm.keys():
                from sklearn.neighbors import NearestNeighbors

                neighs = NearestNeighbors(n_neighbors=n_neighbors + 1)
                neighs.fit(adata.obsm[approx])
                self.indices = neighs.kneighbors_graph(
                    mode="connectivity"
                ).indices.reshape((-1, n_neighbors + 1))
            else:
                from scvelo import Neighbors

                neighs = Neighbors(adata)
                neighs.compute_neighbors(
                    n_neighbors=n_neighbors, use_rep=basis, n_pcs=10
                )
                self.indices = get_indices(
                    dist=neighs.distances, mode_neighbors=mode_neighbors
                )[0]

        self.max_neighs = random_neighbors_at_max

        gkey, gkey_ = f"{vkey}_graph", f"{vkey}_graph_neg"
        self.graph = adata.uns[gkey] if gkey in adata.uns.keys() else []
        self.graph_neg = adata.uns[gkey_] if gkey_ in adata.uns.keys() else []

        if tkey in adata.obs.keys():# 这里修改的内容比较多
            #self.t0 = adata.obs[tkey].copy() 
            self.t0 = adata.obs[tkey].astype("category").copy()
            init = min(self.t0) if isinstance(min(self.t0), int) else 0
            #self.t0.cat.categories = np.arange(init, len(self.t0.cat.categories))
            self.t0 = self.t0.cat.set_categories(np.arange(init, len(self.t0.cat.categories)), rename=True)            
            self.t1 = self.t0.copy()
            #self.t1.cat.categories = self.t0.cat.categories + 1
            self.t1 = self.t1.cat.set_categories(self.t0.cat.categories + 1, rename=True)
        else:
            self.t0 = None

        self.compute_uncertainties = compute_uncertainties
        self.uncertainties = None
        self.self_prob = None
        self.report = report
        self.adata = adata
        
# 原有代码被修改为非并行处理       
    def compute_cosines(self):  
        """TODO."""
        n_obs = self.X.shape[0]

        res = []
        for obs_idx in range(n_obs):
            res.append(self._compute_cosines(obs_idx))

        uncertainties, vals, rows, cols = map(_flatten, zip(*res))

        vals = np.hstack(vals)
        vals[np.isnan(vals)] = 0

        self.graph, self.graph_neg = vals_to_csr(
            vals, rows, cols, shape=(n_obs, n_obs), split_negative=True
        )
        if self.compute_uncertainties:
            uncertainties = np.hstack(uncertainties)
            uncertainties[np.isnan(uncertainties)] = 0
            self.uncertainties = vals_to_csr(
                uncertainties, rows, cols, shape=(n_obs, n_obs), split_negative=False
            )
            self.uncertainties.eliminate_zeros()

        confidence = self.graph.max(1).toarray().flatten()
        self.self_prob = np.clip(np.percentile(confidence, 98) - confidence, 0, 1)

# 原有代码被修改为非并行处理     
    def _compute_cosines(self, obs_idx):  # 删除 queue 参数，适应单线程
    
        vals, rows, cols, uncertainties = [], [], [], []
        
        if self.compute_uncertainties:
            moments = get_moments(self.adata, np.sign(self.V_raw), second_order=True)

        # 将原来的 for obs_id 循环去掉，使用 obs_idx
        
        if self.V[obs_idx].max() != 0 or self.V[obs_idx].min() != 0:
            neighs_idx = get_iterative_indices(
                self.indices, obs_idx, self.n_recurse_neighbors, self.max_neighs
            )

            if self.t0 is not None:
                t0, t1 = self.t0[obs_idx], self.t1[obs_idx]
                if t0 >= 0 and t1 > 0:
                    t1_idx = np.where(self.t0 == t1)[0]
                    if len(t1_idx) > len(neighs_idx):
                        t1_idx = np.random.choice(
                            t1_idx, len(neighs_idx), replace=False
                        )
                    if len(t1_idx) > 0:
                        neighs_idx = np.unique(np.concatenate([neighs_idx, t1_idx]))

            dX = self.X[neighs_idx] - self.X[obs_idx, None]  # 60% of runtime
            if self.sqrt_transform:
                dX = np.sqrt(np.abs(dX)) * np.sign(dX)
            val = cosine_correlation(dX, self.V[obs_idx])  # 40% of runtime

            if self.compute_uncertainties:
                dX /= l2_norm(dX)[:, None]
                uncertainties.extend(
                    np.nansum(dX**2 * moments[obs_idx][None, :], 1)
                )

            vals.extend(val)
            rows.extend(np.ones(len(neighs_idx)) * obs_idx)
            cols.extend(neighs_idx)

        return uncertainties, vals, rows, cols  # 删除 queue 相关处理
    


def _flatten(iterable):
    return [i for it in iterable for i in it]


def velocity_graph(
    data,
    vkey="velocity",
    xkey="Ms",
    tkey=None,
    basis=None,
    n_neighbors=None,
    n_recurse_neighbors=None,
    random_neighbors_at_max=None,
    sqrt_transform=None,
    variance_stabilization=None,
    gene_subset=None,
    compute_uncertainties=None,
    approx=None,
    mode_neighbors="distances",
    copy=False,
    #n_jobs=None,
    #backend="loky",
):

    adata = data.copy() if copy else data
    verify_neighbors(adata)
    if vkey not in adata.layers.keys():
        velocity(adata, vkey=vkey)
    if sqrt_transform is None:
        sqrt_transform = variance_stabilization

    vgraph = VelocityGraph(
        adata,
        vkey=vkey,
        xkey=xkey,
        tkey=tkey,
        basis=basis,
        n_neighbors=n_neighbors,
        approx=approx,
        n_recurse_neighbors=n_recurse_neighbors,
        random_neighbors_at_max=random_neighbors_at_max,
        sqrt_transform=sqrt_transform,
        gene_subset=gene_subset,
        compute_uncertainties=compute_uncertainties,
        report=True,
        mode_neighbors=mode_neighbors,
    )

    if isinstance(basis, str):
        logg.warn(
            f"The velocity graph is computed on {basis} embedding coordinates.\n"
            f"        Consider computing the graph in an unbiased manner \n"
            f"        on full expression space by not specifying basis.\n"
        )


    vgraph.compute_cosines()

    adata.uns[f"{vkey}_graph"] = vgraph.graph
    adata.uns[f"{vkey}_graph_neg"] = vgraph.graph_neg

    if vgraph.uncertainties is not None:
        adata.uns[f"{vkey}_graph_uncertainties"] = vgraph.uncertainties

    adata.obs[f"{vkey}_self_transition"] = vgraph.self_prob

    if f"{vkey}_params" in adata.uns.keys():
        if "embeddings" in adata.uns[f"{vkey}_params"]:
            del adata.uns[f"{vkey}_params"]["embeddings"]
    else:
        adata.uns[f"{vkey}_params"] = {}
    adata.uns[f"{vkey}_params"]["mode_neighbors"] = mode_neighbors
    adata.uns[f"{vkey}_params"]["n_recurse_neighbors"] = vgraph.n_recurse_neighbors

    logg.info("    finished", time=True, end=" " if settings.verbosity > 2 else "\n")
    logg.hint(
        "added \n"
        f"    '{vkey}_graph', sparse matrix with cosine correlations (adata.uns)"
    )

    return adata if copy else None
```

### 6. scvelo RNA速率分析，adipo + seurat UMAP 坐标  
1. python：ob_merge 提取 "PSC_C4", "CTP_C2", "CTP_C3", "VSMC_C2", "Preadipocyte","Adipocyte"  
```
# python
import anndata as ad
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
ob_merge = ad.read('./ob_merge_second.h5ad')
ob_merge.obs['celltype'].value_counts()
scvelo_adipo = ob_merge[ob_merge.obs['celltype'].isin(['PSC_C4', 'CTP_C2', 'CTP_C3', 'VSMC_C2', 'Preadipocyte', 'Adipocyte']), :]
scvelo_adipo.obs['celltype'].value_counts()
```

2. R：seurat 对象提取上述细胞类型重分析后的 UMAP 坐标信息，并转换细胞 id  
```
# R
sub.adipo = readRDS('./sub.adipo.rds') # 经过预处理的 adipo cells
cb <- sub.adipo@meta.data
# umap 坐标
cell_embeddings<-Embeddings(sub.adipo, reduction = "umap")
colnames(cell_embeddings)
cb$UMAP_1 = cell_embeddings[,1]
cb$UMAP_2 = cell_embeddings[,2]
identical(row.names(cb), row.names(cell_embeddings))
scvelo_umap = cb[, c('samplename','UMAP_1', 'UMAP_2')]
# 修改细胞id
scvelo_umap$cell_id<-sapply(row.names(scvelo_umap), 
                        function(x)paste(unlist(strsplit(x, "-"))[1],"x",sep = ""))
# id 头部添加 samplename:
scvelo_umap$cell_id<-paste(scvelo_umap$samplename,
                           scvelo_umap$cell_id,sep = ":")
# 保存
write.csv(scvelo_umap, file = "scvelo_umap.csv", row.names = FALSE)
```

3. seurat UMAP 坐标添加到 anndata 中  
```
# 1. 所需数据
scvelo_umap = pd.read_csv("./scvelo_umap.csv")
scvelo_umap
scvelo_adipo

# 2. 检查cell id是否一致
target = scvelo_umap["cell_id"]
scvelo_adipo = scvelo_adipo[np.isin(scvelo_adipo.obs.index, target)]

# 3. 提取 UMAP 信息，并对其细胞排列顺序
## 二次检查细胞 id
scvelo_adipoe_index = pd.DataFrame(scvelo_adipo.obs.index) # CellID
scvelo_adipoe_index = scvelo_adipoe_index.rename(columns={'CellID': 'cell_id'})
scvelo_umap = scvelo_umap[scvelo_umap['cell_id'].isin(scvelo_adipoe_index['cell_id'])]
## 对齐顺序
scvelo_umap_ordered = scvelo_adipoe_index.merge(scvelo_umap, on='cell_id')
set(scvelo_umap_ordered[['cell_id']]) == set(scvelo_adipoe_index[['cell_id']])  # 检查顺序
# 仅保留 UMAP_1 和 UMAP_2 两列
scvelo_umap_ordered = scvelo_umap_ordered[['UMAP_1', 'UMAP_2']]

# 4. 添加 UMAP 信息
scvelo_adipo.obsm['X_umap'] = scvelo_umap_ordered.values
scvelo_adipo.obsm
```

4. RNA 速率分析
**最好先跑一遍示例数据，看看能不能跑通**


adata = scv.datasets.pancreas()
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
# 
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120)
plt.savefig('embedding.png', dpi=300, bbox_inches='tight')
plt.close()

scv.pl.velocity_embedding_stream(adata, basis='umap')
plt.savefig('embedding_stream.png', dpi=300, bbox_inches='tight')
plt.close()





adata = ad.read('./example.h5ad')
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120)
plt.savefig('embedding_example.png', dpi=300, bbox_inches='tight')
plt.close()

scv.pl.n

(adata, basis='umap')
plt.savefig('embedding_stream_example.png', dpi=300, bbox_inches='tight')
plt.close()





scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120)

ob_merge = ad.read('./ob_merge_second.h5ad')

# 默认图片参数
scv.set_figure_params()
scv.logging.print_version()
# scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

# 提取子集
ob_merge.obs['celltype'].value_counts()
example = ob_merge[ob_merge.obs['celltype'].isin(['PSC_C3', 'VSMC_C1', 'VSMC_C2']), :]
del ob_merge
ob_merge

## scvelo RNA速率分析
scv.pp.filter_and_normalize(example)
scv.pp.moments(example)
scv.tl.velocity(example, mode = "stochastic")
scv.tl.velocity_graph(example)



adata = scv.datasets.pancreas()
example = scv.datasets.pancreas()
adata
adata = example

scv.pl.proportions(adata)
plt.savefig('proportions_plot_example.png', dpi=300, bbox_inches='tight')
plt.close()


adata.write('./example.h5ad')



adata = ad.read('./example.h5ad')
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata)
scv.tl.velocity_graph(example)


















# spliced/unspliced的比例
scv.pl.proportions(adata)

# plot
scv.pl.velocity_embedding(adata, basis = 'umap')
scv.pl.velocity_embedding_stream(adata, basis = 'umap')

## scvelo RNA速率分析
# 预处理1
scv.pp.filter_and_normalize(ob_merge)
scv.pp.moments(ob_merge)
scv.tl.velocity(ob_merge, mode = "stochastic")
scv.tl.velocity_graph(ob_merge)

# 预处理2
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)



# spliced/unspliced的比例
scv.pl.proportions(adata)

# 可视化
scv.pl.velocity_embedding(adata, basis = 'umap')


scv.pl.velocity_embedding_stream(adata, basis = 'umap')









