# Pan-TAG v0.1

Pan-TAG 是一个基于 Snakemake 的泛基因组串联重复基因簇（tandemly arrayed genes，TAG）分析工具。它从同一物种的多套基因组组装和注释出发，完成单样本 TAG 识别、可选假基因召回、跨样本合并和全样本 refinement；还可按需提取序列、进行 Pan-cluster 内序列聚类，并绘制 PCoA 与基因排列结构图。

核心输出：

```text
results/local_tags/<sample_id>_TAGs.tsv
results/global/pangenome_TAGs.tsv
results/global/refined_TAG_count_matrix.tsv
```

当前版本支持 Linux x86-64，软件环境由 Pixi 管理。

## 1. 下载与安装

```bash
cd /path/to/software
git clone https://github.com/rongjin625/Pan-TAG.git
cd Pan-TAG
pixi install
```

如果要让 Pan-TAG 为非参考样本运行 Liftoff，再安装一次独立环境：

```bash
pixi install --environment liftoff
```

环境安装应在可联网节点完成，不要在离线计算节点重复执行。检查安装：

```bash
pixi run pantag --version
pixi run snakemake --version
pixi run diamond version
pixi run gffread --version
```

定制版 MCScanX 无需单独安装，第一次运行工作流时会从项目内源码自动编译。

## 2. 建立数据集运行目录

代码和环境可以由多个数据集共用；每个数据集应使用独立运行目录：

```bash
CODE_ROOT=/path/to/software/Pan-TAG
RUN_ROOT=/path/to/pantag_runs/my_dataset

mkdir -p "$RUN_ROOT/config"
cd "$RUN_ROOT"
```

后续命令都应从 `$RUN_ROOT` 执行。输入 FASTA/GFF 建议填写绝对路径。

## 3. 准备 `samples.tsv`

创建 `$RUN_ROOT/config/samples.tsv`。它必须是制表符分隔文件，表头和列顺序固定为：

```text
sample_id	fasta	gff
```

共同要求：

- 第一条数据必须是参考样本；
- 参考样本必须同时提供 FASTA 和 GFF；
- `sample_id` 必须唯一；
- FASTA 序列名必须与相应 GFF 第一列 `seqid` 一致；
- 一次分析只放同一物种的样本。

### 使用已有注释

所有样本都填写 GFF：

```tsv
sample_id	fasta	gff
reference	/data/reference.fa	/data/reference.gff3
sample01	/data/sample01.fa	/data/sample01.gff3
sample02	/data/sample02.fa	/data/sample02.gff3
```

此时使用 `annotation.mode: provided`。

### 由 Pan-TAG 运行 Liftoff

只有首行参考样本填写 GFF，其他样本的第三列留空：

| sample_id | fasta | gff |
|---|---|---|
| reference | `/data/reference.fa` | `/data/reference.gff3` |
| sample01 | `/data/sample01.fa` | 留空 |
| sample02 | `/data/sample02.fa` | 留空 |

此时使用 `annotation.mode: liftoff`。Pan-TAG 会先用参考 GFF 建立数据库，再对每个非参考样本运行 Liftoff。一次运行不能混用两种注释模式。

Liftoff 的 `-f` 类型列表固定读取 `$CODE_ROOT/config/liftoff_features.txt`；模板当前包含 `gene` 和 `pseudogene`。只有确实需要转移其他顶层 feature type 时才修改这个项目文件。

## 4. 生成 `config.yaml`

不要直接使用仓库模板中的示例数据路径。根据 `samples.tsv` 生成本次运行的配置。

使用已有 GFF：

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag configure \
  --samples "$RUN_ROOT/config/samples.tsv" \
  --annotation-mode provided \
  --pseudogene-enabled \
  --template "$CODE_ROOT/config/config.yaml" \
  --output "$RUN_ROOT/config/config.yaml" \
  --results-dir "$RUN_ROOT/results" \
  --logs-dir "$RUN_ROOT/logs" \
  --work-dir "$RUN_ROOT/work"
```

由 Pan-TAG 运行 Liftoff：

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag configure \
  --samples "$RUN_ROOT/config/samples.tsv" \
  --annotation-mode liftoff \
  --pseudogene-enabled \
  --template "$CODE_ROOT/config/config.yaml" \
  --output "$RUN_ROOT/config/config.yaml" \
  --results-dir "$RUN_ROOT/results" \
  --logs-dir "$RUN_ROOT/logs" \
  --work-dir "$RUN_ROOT/work"
```

不需要假基因召回时，把 `--pseudogene-enabled` 换成 `--no-pseudogene-enabled`。查看全部可调选项：

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag configure -h
```

### 新物种的 GFF 设置

默认模板兼容常见 GENCODE、Ensembl 和标准 GFF3：

```yaml
annotation:
  gene_feature_types: [gene]
  pseudogene_feature_types: [gene, pseudogene]
  transcript_feature_types: [mRNA, transcript]
  gene_id_keys: [gene_id, ID]
  gene_name_keys: [gene_name, Name, gene_symbol, gene, geneName]
  transcript_id_keys: [transcript_id, ID]
  transcript_parent_gene_keys: [gene_id, Parent]
  biotype_keys: [gene_type, gene_biotype, biotype]
  coding_filter_mode: auto
  coding_values: [protein_coding]
  excluded_attribute_terms: [readthrough_gene, readthrough]
  pseudogene_values: [pseudogene, processed_pseudogene, unprocessed_pseudogene]
  strip_id_prefixes: ["gene:", "transcript:"]
  fallback_gene_name_to_id: true
```

首次处理新注释来源时，应检查 GFF 第三列和第九列的 ID、Name、Parent 与 biotype 写法，并相应修改生成的 `config.yaml`。假基因只有在 feature type 属于 `pseudogene_feature_types` 且 biotype 匹配 `pseudogene_values` 时才会被读取；默认不读取 `pseudogenic_transcript`。

`coding_filter_mode` 可设为：

- `auto`：有可识别的 coding biotype 时按 biotype 筛选，否则按有效翻译蛋白筛选；
- `biotype`：强制按 `coding_values` 筛选；
- `translated_protein`：只根据能否获得有效翻译蛋白筛选。

## 5. Dry-run

先检查核心 DAG：

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" snakemake all \
  --directory "$RUN_ROOT" \
  --snakefile "$CODE_ROOT/workflow/Snakefile" \
  --configfile "$RUN_ROOT/config/config.yaml" \
  --cores 20 \
  --dry-run \
  --printshellcmds
```

`provided` 模式不应出现 `liftoff`；`liftoff` 模式应为每个非参考样本生成一个 Liftoff job。

## 6. 运行核心分析

### 只运行核心流程

运行到最终 Pan-cluster 表和 refinement 计数矩阵，不提取序列、不聚类：

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" snakemake all \
  --directory "$RUN_ROOT" \
  --snakefile "$CODE_ROOT/workflow/Snakefile" \
  --configfile "$RUN_ROOT/config/config.yaml" \
  --cores 20 \
  --rerun-incomplete \
  --printshellcmds
```

### 运行下游序列分析

核心流程完成后，以下命令会按需补齐所需的核心结果和序列，并对所有 Pan-cluster 进行序列聚类、PCoA、结构图和汇总：

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag cluster \
  --config "$RUN_ROOT/config/config.yaml" \
  --all \
  --cores 20 \
  --printshellcmds
```

提交到调度系统时，`--cores` 不应超过任务实际申请的 CPU 数量。

## 7. 按需提取与聚类

核心流程完成后，可以随时选择感兴趣的 Pan-cluster，不需要修改 `config.yaml`。

只提取序列：

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag extract \
  --config "$RUN_ROOT/config/config.yaml" \
  --tag Pan_cluster_0132 \
  --tag Pan_cluster_0433 \
  --cores 8
```

聚类并绘图：

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" pantag cluster \
  --config "$RUN_ROOT/config/config.yaml" \
  --tag Pan_cluster_0132 \
  --tag Pan_cluster_0433 \
  --cores 8
```

将重复的 `--tag` 换成 `--all` 即处理全部 Pan-cluster。两个命令都支持先加 `--dry-run --printshellcmds` 查看计划。

`pantag cluster` 会自动补齐缺少的核心结果和序列。对同一次选择，所有簇先完成序列聚类；随后仅为 `status=ok` 且具有有效 PCoA 坐标的簇绘制结构图。序列不足或没有可用 k-mer 的簇会保留状态记录，但不会生成有效的结构图。

## 8. 主要参数

参数都保存在生成的 `config.yaml` 中：

| 配置 | 含义 | 模板值 |
|---|---|---:|
| `liftoff.threads` | 每个 Liftoff job 的线程数 | 8 |
| `tag.diamond_threads` | 每个 DIAMOND blastp job 的线程数 | 4 |
| `tag.diamond_max_target_seqs` | 每个 query 最多报告的非自身候选数 | 5 |
| `tag.blastp_*` | blastp identity/query/subject coverage 阈值 | 20 |
| `pseudogene.search_flank_bp` | 假基因搜索侧翼范围 | 50000 |
| `pseudogene.blastx_*` | blastx identity/query/subject coverage 阈值 | 20 |
| `flanks.gene_count` | TAG 每侧参与 QC/共线性的基因数 | 10 |
| `merge.member_overlap_threshold` | 跨样本成员重叠阈值 | 0.3 |
| `merge.require_same_seqid` | 是否限制候选 TAG 位于同名 seqid | false |
| `refine.flank_overlap_threshold` | refine 侧翼匹配阈值 | 0.3 |
| `downstream.cluster.k` | canonical k-mer 长度 | 31 |
| `downstream.cluster.ubiquitous_threshold` | 在序列集合中出现比例达到该值的 k-mer 将被过滤 | 0.95 |
| `downstream.cluster.max_clusters` | 最大分组数；0 表示自动 | 0 |
| `downstream.cluster.silhouette_tolerance` | 优先较少分组时允许的 silhouette 差值 | 0.05 |
| `downstream.cluster.min_sequences` | 启动聚类的最少序列数 | 10 |
| `downstream.cluster.pcoa` | 是否输出 PCoA 图 | true |

DIAMOND blastp 使用 query coverage 和 subject coverage，并通过 `--no-self-hits` 排除自身命中。单样本 MCScanX 输入会跳过完全包含于前一保留 gene 区间内的嵌套基因，数量记录在 `<sample>_prepare_qc.tsv` 的 `nested_genes_skipped` 中。

## 9. 输出文件

```text
results/
├── local_tags/
│   └── <sample_id>_TAGs.tsv
├── global/
│   ├── pangenome_TAGs.tsv
│   └── refined_TAG_count_matrix.tsv
└── downstream/
    ├── global_cluster_summary.tsv
    └── Pan_cluster_XXXX/
        ├── sequences/
        │   └── Pan_cluster_XXXX_sequences.fa
        └── clusters/
            ├── Pan_cluster_XXXX_pcoa.pdf
            ├── Pan_cluster_XXXX_pcoa_coordinates.tsv
            ├── Pan_cluster_XXXX_TAG_structure.pdf
            └── Pan_cluster_XXXX_TAG_structure_collapsed.pdf
```

- `<sample_id>_TAGs.tsv`：单样本 cleaned TAG，包含成员、坐标、密度和两侧基因；存在召回假基因时增加假基因列。
- `pangenome_TAGs.tsv`：最终 Pan-cluster 总表，包含成员名称/ID、参考 seqid、合并后的左右侧翼基因以及 refinement 所需的成员信息。
- `refined_TAG_count_matrix.tsv`：Pan-cluster × sample 计数矩阵及检测样本数、无 gap 样本数、均值、标准差和变异系数。
- `Pan_cluster_XXXX_sequences.fa`：该 Pan-cluster 中无 gap、非 span 记录的序列，标题同时包含 sample、gene name 和 gene ID。
- `Pan_cluster_XXXX_pcoa.pdf`：相同坐标的双面板图；左图按原始 gene name 着色，右图按序列分组 `Pan_cluster_XXXX-SGXX` 着色。
- `Pan_cluster_XXXX_pcoa_coordinates.tsv`：逐序列的 PCoA 坐标、gene 信息、Sequence Group 和固定颜色。
- `Pan_cluster_XXXX_TAG_structure.pdf`：按样本展示 TAG 基因排列；目标基因使用 Sequence Group 颜色，侧翼基因显示为灰色。仅对序列聚类成功且具有有效 PCoA 坐标的 Pan-cluster 生成。
- `Pan_cluster_XXXX_TAG_structure_collapsed.pdf`：将具有相同 TAG 组成的样本折叠后展示，并标注对应样本数量。
- `global_cluster_summary.tsv`：所有已完成聚类的状态和主要指标汇总。

计数矩阵的样本单元格：

| 值 | 含义 |
|---|---|
| 正整数 | 唯一有效区块，且目标区间不含 N |
| `N*` | 唯一有效区块，但目标区间含 N |
| `N#` | 存在多个侧翼支持区块（span） |
| `NA` | 未找到目标基因或没有通过侧翼 QC 的区块 |

仅当一个候选簇的 span 样本比例严格大于 0.9 时，refine 才将该候选簇从最终结果删除。假基因计入样本拷贝数和汇总统计。

当 `outputs.debug: true` 时，额外发布：

```text
results/global/local_to_global.tsv
```

日志写入 `logs_dir`，中间文件与断点续跑状态写入 `work_dir`。首次使用新注释时重点检查：

```text
work/<sample>/<sample>_prepare_qc.tsv
work/<sample>/<sample>.map
work/downstream/TAG_sequence_manifest.tsv
```

`TAG_sequence_manifest.tsv` 是下游序列提取和聚类使用的内部清单，记录通过 refinement 且可提取的 TAG 成员坐标。

## 10. 断点续跑与多数据集并行

任务失败后直接用原命令重新运行。Snakemake 会在 `--rerun-incomplete` 下重建不完整文件，并复用仍然有效的结果。

若异常退出后遗留锁，必须先确认该运行目录没有其他 Snakemake 任务，再执行：

```bash
pixi run --manifest-path "$CODE_ROOT/pixi.toml" snakemake \
  --directory "$RUN_ROOT" \
  --snakefile "$CODE_ROOT/workflow/Snakefile" \
  --configfile "$RUN_ROOT/config/config.yaml" \
  --unlock
```

不同数据集可以共用 `CODE_ROOT` 和 Pixi 环境，但必须使用不同的 `RUN_ROOT`、`results_dir`、`logs_dir` 和 `work_dir`。不要让两个任务共享 `work_dir`。

常见提示：

- `Nothing to be done`：目标已经存在且为最新状态，不是失败；
- `IncompleteFilesException`：使用带 `--rerun-incomplete` 的原命令重跑，不要随意把未验证文件标记为完整；
- `Directory cannot be locked`：仍有任务运行或上次异常退出留下锁；
- 修改 `config.yaml` 后任务数增加：先执行 `--dry-run`，Snakemake 可能因输入或参数变化重跑受影响步骤。

## 11. 开发测试

```bash
cd "$CODE_ROOT"
PYTHONPATH=src pixi run python -m unittest discover -s tests -v
```

## Citation

Pan-TAG v0.1 仍在开发中，论文与正式引用信息将在发布后补充。

## License

当前仓库尚未添加项目级 `LICENSE`。定制 MCScanX 的上游许可见 `src/pantag/mcscanx/LICENSE`。
