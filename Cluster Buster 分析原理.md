## [Cluster Buster](https://bu.wenglab.org/cluster-buster/index.html)
该方法用于识别某段DNA序列可能存在的motif簇（motif cluster）和可能性，及其具体富集到的motif基序和得分


### 分析背景：
1. 转录的增强或抑制是由 **多个转录因子结合位点的集合（clusters of transcription factor binding sites）** 调控的，也就是说这一段区域可能同时结合多个转录因子，具有多个转录因子的motif序列，因此称为motif cluster
2. 目前有两种方法评估motif基序与DNA序列关系：一种是直接对某一特定DNA序列窗口出现的motif基序计数，这种方法更加只管；另一种是利用**概率模型算法**（不懂），也就是 cbust 软件使用的方法
3. 两个输入文件：FASTA格式目标DNA序列和motifs集合
4. 可调节参数
5. 输出数据
   * motif cluster在某段DNA序列上的分布和得分
   * 单个motif cluster 所包含的全部motif基序和得分


### 分析原理
**针对给定的某一DNA片段或多个目标DNA片段**
1. 执行前向算法
   * 前向算法（Forward algorithm）：从第一个核苷酸开始，连续计算到第 i 个核苷酸的每个子序列的对数似然得分（log likelihood score s[i]）。这个得分帮助识别序列中潜在的功能区域。
   * 追踪得分增加最大的子序列：标记那些得分增加s[b]−s[a] 最大的子序列(a,b)，这里 a 是子序列的开始点，b 是结束点。算法确保这些子序列不会与其他得分更高的子序列重叠，从而找到局部最优的区域。
2. 使用后向算法（backward algorithm）细化起始点
   * 每个识别的子序列 (a,b)中，b 点也就是序列终点被认为是可靠的，而起始点 a 可能不准确。
   * 通过后向算法：从 b 点开始，向 a 反向执行后向算法，直到接近 a 点前。这一步骤旨在优化和精确确定最佳的起始点a，进一步提高识别的精度。
3. 移除重叠的高得分子序列
   * 使用贪婪算法（greedy algorithm）移除那些与得分更高的子序列重叠的区间。这意味着如果有两个或更多重叠的子序列，算法将保留得分最高的那个，删除其他得分较低的子序列。
4. **通过这种方法，Cluster-Buster 能够有效地从复杂的基因组数据中筛选出重要的生物学信号，这对于理解基因调控机制和发现新的调控元件至关重要**

### 输入文件：FASTA 文件
```
>FirstSeq
cgcccctacaagtatcaggtagagtgctccacaggtgaggtaaaacgaacgctccatgggggcgcctatcgcattcaattaaaatacgacacattacgg
atggtggttcccgggacatgggtaacgtcagtagcacaccaccagcgtttatgttgcacttacgggactaagttacttaaactgttcaggagatacccc
>2ndSeq
gaaaagtaaaatcgtcccaacacattcgcacgaccttcaggtggatggcgaccgctcctaggaccggtgtgcggcgtagccccgggttgacaaaaagga
caccttcaaagcgtagaagatgccaatggtctgtgcatttaagggacctgataaccgcacccatagtgacctataaccccgaatgttaccatctactca
>3rdSeq
aataagcgaagagcttaaagtcatcgggtcagagatgatcagtcaaggatgcatcacctcatggaagcctaacgtgagtgcgggccctctatcagcccg
atgagatgggagcactcgtgttctagatgtttaccgtgtttcgaagatttattaaccaggaggttttcgctaaccgagtcctaatgcgatctatctaca
......
>nthSeq
tctctcttagtagattttaacagcctgatatggcggcttcacaattggtgacggcttccaggtgtccacgtcctgacggtcttgaaatctttccaatgc
```
### 输入文件：motif 文件；位置权重矩阵；PWM
从5'-3'依次为A, C, G, T，大小写均可，其他无法识别字符被转为n并跳过
```
>element1
0  4 2 14
12 0 0 8
8  0 1 11
20 0 0 0
>element2
13 1 1 5
```
### 输出文件：motif cluster name, location, strand, score, and sequence 及其所包含的motif信息
```
CLUSTER 16
>gi|312167|emb|X70518.1|RNPROL R. norvegicus prolactin gene (5' region)  (2014 bases)
Location: 811 to 1045
Score: 6.06
CCAGGTCATCTGTCagtccaaattcagaaacagtaaagccaaaactaaaggtCACAAGCTGCTTcagatgaatgaatccc
caaattaaagaaagtcatcagcaacttcattattattcaccataatgacatcatttaggaaatctctaaaacatgagtgg
aactttggagtgcattaaaaaatgcatttTTGTCACTATGTCCTagagtgctttggGGTCAGAAGAGGCAGGCAG
V$SF1_Q6                 811 -  818   -   1.44       tgacctgg
V$E12_Q6                 814 -  824   -   3.36       gacagatgacc
Myf                      863 -  874   -   1.68       aagcagcttgtg
ERE                     1000 - 1014   -   0.729      aggacatagtgacaa
V$TANTIGEN_B            1027 - 1045   +   2.3        ggtcagaagaggcaggcag
V$CACBINDINGPROTEIN_Q6  1035 - 1043   +   1.52       gaggcaggc
```




