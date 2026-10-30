# 随Git分发的公开图谱

这4个图谱是公开研究资源，不是LS/C3/GB患者影像。图像保持原始发布字节；TSV由原XML转换，编号规则见PUBLIC_MANIFEST.json和docs/ATLAS_RESOURCES.md。原XML随仓库保存在resources/atlas_sources。

| 文件组 | 作者/提供者及引用 | 许可 |
|---|---|---|
| HarvardOxford皮层及皮层下 | Harvard Center for Morphometric Analysis及FSL；Desikan等，NeuroImage 2006，31:968–980；其他引用见FSL图谱页 | CC BY-SA 4.0 |
| Cerebellum MNIflirt | Jörn Diedrichsen、Narender Ramnani及合作者；Diedrichsen等，NeuroImage 2009，46:39–46 | CC BY-SA 4.0 |
| JHU20概率束 | Susumu Mori及Johns Hopkins University；Hua等，NeuroImage 2008，39:336–347 | FSL非商业许可 |

发布来源及许可于2026-10-30核对：

- https://fsl.fmrib.ox.ac.uk/fsl/docs/other/datasets.html
- https://fsl.fmrib.ox.ac.uk/fsl/docs/license.html
- https://creativecommons.org/licenses/by-sa/4.0/
- 侧化HO原始分发包：https://www.nitrc.org/frs/download.php/9902/HarvardOxford.tgz

完整FSL主许可保存在resources/licenses/FSL-LICENSE.txt；各图谱适用的许可范围由上述FSL官方许可页面明确列出。CC BY-SA适用于HO和小脑及其派生标签，不以本项目许可覆盖；JHU继续遵守FSL非商业条件。本仓库免费分发数据及原XML，没有复制或修改FSL编译程序。

在仓库根目录校验：

```bash
sha256sum -c resources/atlases/SHA256SUMS
```

主配置已固定4个NIfTI的SHA256。正常运行无需重装FSL或重新下载图谱。atlas_manifest.json是本机准备时的完整路径记录，不进Git；PUBLIC_MANIFEST.json保留公开来源、相对路径、哈希、版本及QC。
