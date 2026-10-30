# 急性病灶图谱资源

资源准备与患者分析分开。脚本只读取已安装FSL、获取缺少的公开侧化图谱，并生成本项目所需的4个2 mm图谱及标签表；不处理患者影像，不修改源FSL文件。

```bash
.venv/bin/python scripts/prepare_atlases.py \
  --config config/project.yaml \
  --sources config/atlas_sources.yaml
```

主config确定图谱输出位置、模板、概率尺度和标签排除规则；`atlas_sources.yaml`确定安装版FSL位置、版本、源文件、XML索引规则、下载策略、缓存和来源清单位置。默认拒绝覆盖不同内容的现有资源。已相同的资源保持原内容及修改时间。

## 本机已准备的资源

2026-10-30完成实际文件准备和NiBabel检查：

| 图谱 | 实际来源 | 使用的区域数 | 编号规则 |
|---|---|---:|---|
| Harvard–Oxford侧化皮层 | Nilearn 0.13.1获取的FSL公开预侧化版本，原文件逐字节复制 | 96 | 输出NIfTI整数标签对应TSV index；0为背景 |
| Harvard–Oxford皮层下及脑干 | 已安装FSL 6.0.7.18，逐字节复制 | 15 | XML index＋1；按主配置排除泛灰质、白质、脑室和背景 |
| JHU20概率纤维束 | 已安装FSL 6.0.7.18，逐字节复制 | 20 | XML index直接对应从0开始的4D体积索引 |
| Diedrichsen MNIflirt小脑 | 已安装FSL 6.0.7.18，逐字节复制 | 28 | XML index＋1；0为背景 |

合计159个区域/纤维束特征。它们不是互斥的完整脑分区；概率纤维束之间可以重叠，部分病灶可不被图谱覆盖。每个患者仍需输出未覆盖病灶比例。

4个图谱空间尺寸均为91×109×91，JHU额外有20个概率体积；体素大小均为2 mm，NIfTI空间单位均明确为mm。模板声明为`MNI152NLin6Asym`。这一声明来自资源定义，不是通过“都是MNI、都是2 mm”推断出来的患者配准证明。

## 侧化方法和已发表边界

目前实际使用公开预侧化图谱，未执行自定义左右半球划分。Nilearn现行`fetch_atlas_harvard_oxford(..., symmetric_split=True)`会请求预侧化文件，并非一定在本地非侧化图谱上计算划分；脚本因此优先明确请求`cortl-maxprob-thr25-2mm`，下载缓存位于项目临时目录，不写入原FSL安装。

若公开资源获取失败且配置允许，备用方案才调用已锁定Nilearn自身的`_compute_symmetric_split`，并记录软件版本、辅助函数源代码SHA256及下载失败原因。备用输出会明确标为派生图谱，不冒充公开预侧化原文件。

侧别核查基于MNI世界坐标质心和体素比例，不删改已发表图谱的跨中线边界。实际96个侧化标签中，只有`Left Occipital Pole`包含少量对侧体素：2746个体素中20个处于x=+2 mm，比例0.7283%，质心x=-16.311 mm。其余95个标签没有对侧体素。配置要求质心位于正确半球、对侧比例不超过5%，逐标签数值均写入来源清单。

## 完整性与传递

`resources/atlases/atlas_manifest.json`记录源文件及输出文件SHA256、FSL/Nilearn/NiBabel版本、标签编号规则、空间单位、仿射矩阵、许可和侧别核查结果。三个本地FSL图谱和公开侧化皮层图谱的输出均与各自源图像字节一致。源FSL图像与XML在准备前后重新计算哈希，确认未被修改。

运行前严格拒绝未知、米或微米单位，不悄悄重解释头信息。缺少资源或哈希不一致时明确报告，不能用其他图谱临时替代。

2026-10-30 Git交付调整后，4个公开图谱、标签、原XML和许可已加入Git明确名单。目标工作站无需为运行本项目重新安装整套FSL。主config使用相对路径并固定NIfTI SHA256；`resources/atlases/SHA256SUMS`校验随仓库分发的文件。

本机完整准备清单`atlas_manifest.json`仍保留在私有备份；公开版`PUBLIC_MANIFEST.json`保留来源、相对路径、哈希和QC，去除本机根路径。图像和标签的科学内容未变化。

## 来源与许可

- [Nilearn Harvard–Oxford接口](https://nilearn.github.io/stable/modules/generated/nilearn.datasets.fetch_atlas_harvard_oxford.html)
- [FSL图谱定义与引用](https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html)
- [FSL许可说明](https://fsl.fmrib.ox.ac.uk/fsl/docs/license.html)

主配置将Harvard–Oxford和小脑图谱标为CC BY-SA 4.0、JHU标为FSL非商业许可。图谱许可与本项目代码许可分别记录，物理资源包保留来源信息；资源准备不意味着改变原许可。
