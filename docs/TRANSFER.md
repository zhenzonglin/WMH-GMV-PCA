# Git与物理传递清单

更新：2026-10-30。仓库：https://github.com/zhenzonglin/WMH-GMV-PCA 。本次最大化随Git分发项目可复用内容；旧版code/atlases分离zip保留为历史交付，新克隆已经包含公开图谱。

## 随Git传递

| 内容 | 路径 |
|---|---|
| Python/R分析代码、调度及管理入口 | src/、R/、scripts/、workflow/ |
| 可移植配置与严格schema | config/，不含config/local |
| 统计方案、字段字典、操作和验收记录 | docs/、README.md |
| 完整依赖版本锁及安装脚本 | requirements.lock、renv.lock、scripts/setup* |
| 4个公开图谱、4个标签表、原XML、许可及哈希 | resources/ |
| 测试和合成数据生成方法 | tests/、scripts/make_synthetic.py |
| 合成验收摘要及截图 | docs/verification/synthetic_v2/ |

公开图谱原NIfTI与标签合计约270 KiB；含完整本机来源清单时约312 KiB，远低于GitHub单文件限制。其内容不变，只调整传递方式。JHU继续保留非商业使用限制；HO和小脑继续采用CC BY-SA 4.0，详见resources/atlases/README.md。

普通Git对超过50 MiB的单文件警告，超过100 MiB的单文件拒绝；Git LFS或Releases是其他分发途径。本项目本次纳入Git的文件无需LFS。依据：https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github 。目录体积大并不等于其中每个文件不能上传，应区分技术限制和是否适合纳入版本库。

## 不进入Git及物理备份

| 内容 | 原因 | 本次处理 |
|---|---|---|
| 真实患者影像、表格及患者模型 | 不适合进入公开项目仓库 | 本地没有真实队列输入，不能伪造或打包；由远程私有存储另行传递 |
| 完整合成运行、失败历史、日志、输入、私有路径配置 | 保留研究记录，避免重复输出充满代码历史 | local-records压缩包 |
| .venv、.Rlib安装目录 | 约894 MiB，依赖平台和原绝对路径；锁文件可重建 | environment-snapshot压缩包，作为本机备份 |
| 临时计算诊断 | 属于验证记录 | 并入local-records，排除可重新下载的图谱缓存 |
| .git、缓存、旧dist压缩包、编译缓存 | 可再生成或重复备份 | 不重复打包 |

物理包：

- `WMH-GMV-PCA-local-records-20261030.tar.gz`：项目目录内的data、runs、logs、state、config/local、.verification与诊断tmp；只包含本机现有资料，当前均为软件合成测试相关内容。还保存原始atlas_manifest.json与旧名册路径参考。
- `WMH-GMV-PCA-environment-snapshot-20261030.tar.gz`：本机.venv和.Rlib、环境指纹及锁文件。不包含所有外部系统依赖；Python解释器符号链接和R的其他库路径在说明中列出。这不是可保证跨机器直接启动的离线安装器。

两个包为可选备份。正常远程部署需要Git克隆、按锁文件安装依赖、配置真实队列路径，不需要把本机历史运行或安装目录覆盖进新项目。

## 接收与核验

先将物理包和SHA256SUMS.txt放在同一目录：

```bash
sha256sum -c SHA256SUMS.txt
mkdir -p history_copy
tar -xzf WMH-GMV-PCA-local-records-20261030.tar.gz -C history_copy
```

每个包有独立顶层目录、MANIFEST.json与FILE_MANIFEST.sha256。进入解压后的顶层目录执行`sha256sum -c FILE_MANIFEST.sha256`核对所有普通文件；MANIFEST.json另外记录符号链接目标。不要对未检查的旧环境执行activate，也不要将旧主机的PID/锁记录当成当前状态。

首次克隆：

```bash
git clone https://github.com/zhenzonglin/WMH-GMV-PCA.git
cd WMH-GMV-PCA
sha256sum -c resources/atlases/SHA256SUMS
bash scripts/setup_project.sh --config config/project.yaml
# 设置实际paths和临床字段映射后：
.venv/bin/python scripts/audit_inputs.py --config config/project.yaml
```

真实数据或私有配置应按当地资料管理要求传递。原工作站路径参考位于local-records包的`config/local/workstation-original-20261030/`，它仍需远程核实，不能代替输入审计。
