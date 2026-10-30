# 第三方软件与图谱

本项目调用成熟软件，不复制或修改其统计求解器。完整组合流程仍需自身校准、稳定性及外部验证。

| 用途 | 上游 |
|---|---|
| GAMLSS与分布 | https://github.com/gamlss-dev/gamlss ; https://github.com/gamlss-dev/gamlss.dist |
| PCA及投射 | https://github.com/husson/FactoMineR |
| PCA维数交叉验证 | https://github.com/husson/missMDA |
| 有序回归 | https://github.com/runehaubo/ordinal |
| 影像读取及图谱操作 | https://github.com/nipy/nibabel ; https://github.com/nilearn/nilearn |
| 调度 | https://github.com/snakemake/snakemake |
| 环境 | https://github.com/rstudio/renv ; https://github.com/astral-sh/uv |

Snakemake安全停止语义按验证的9.27.1版本：仅向调度器发送SIGTERM进入排空；不采用Ctrl-C取消语义。

图谱说明：https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html 。许可：https://fsl.fmrib.ox.ac.uk/fsl/docs/license.html 。HO及小脑图谱按其官方CC BY-SA条款，JHU按官方FSL非商业条款核查。软件许可证不自动覆盖图谱数据。

本仓库附4个公开图谱及原XML，按原许可免费分发，不附患者影像。署名、引用和许可范围见`resources/atlases/README.md`；完整FSL许可见`resources/licenses/FSL-LICENSE.txt`。版本、来源、标签、许可与哈希同时写入运行记录。GM119标准空间显示图谱未取得前，使用载荷表/热图，不生成冒充精确配准的脑图。
