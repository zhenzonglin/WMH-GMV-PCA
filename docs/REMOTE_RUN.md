# 远程工作站运行说明

## 目录与环境

将代码放到工作站的独立项目目录。更新`config/project.yaml`的`project_root`以及`config/paths.yaml`。每台工作站的私有配置建议放在被Git忽略的`config/local/<machine>/`，完整复制5份配置后调整；主文件仍按相同规则引用其余4份。

普通克隆位于任意目录时，默认`project_root: ..`无需修改。主配置若复制到`config/local/<machine>/`，可设`project_root: ../../..`或项目绝对路径。4个公开图谱已经包含于仓库，先在根目录执行`sha256sum -c resources/atlases/SHA256SUMS`；不再需要单独导入旧版图谱包。

已测试Python 3.11、R 4.5.2、Snakemake 9.27.1。源码安装R包需要C/C++/Fortran编译器、make以及CMake；CMake包含于已验证Python锁文件。远程R环境不同时应先重建对应版本，不能将未经验证的运行写成已验证结果。

入口安装：`bash scripts/setup_project.sh --config config/project.yaml`。需要uv在PATH、配置的Rscript存在。正式依赖依据锁文件恢复；安装日志为配置`log_root`下的`install_r.log`。

## 输入核对

运行`.venv/bin/python scripts/audit_inputs.py --config config/project.yaml`。

必需影像数据为三队列名册、年龄/性别、WMH20、GM119及真实ICV。可配置显式TSV；也支持`raw_records.jsonl.gz`中`features.wmh.raw_ml`、`features.t1.gm119_ml`及`dlicv_icv.icv_label702_ml`。原始GM标签保留，不使用旧macro20替代。已有分区表存在多个候选时必须指定路径。

WMH全脑总量不足时，仅从已有校正二值掩膜统计；要求明确mm单位，不重新处理掩膜，不把分区之和默认为全脑量。临床模型缺总量时报告缺项。

病灶掩膜需要`lesion_template`声明，与图谱模板一致。统一为MNI152NLin6Asym的声明必须有上游配准依据；相同尺寸或文件名不能代替此依据。只允许同模板的网格重采样。掩膜及图谱要求NIfTI空间单位为mm。原生病灶需要上游提供已配准版本。

图谱路径和标签表在`paths.yaml`配置；标签表至少包含`index,name`，概率束索引是0开始的4D volume编号。HO皮层下排除通用WM、GM和脑室标签；不能把这些总组织标签当局部灰质区域。JHU20概率图不是JHU48标签图。

## 启动和监测

在项目目录执行：

```bash
bash scripts/start_project.sh --config config/project.yaml
.venv/bin/python scripts/monitor_project.py --config config/project.yaml
```

脚本打印run_id、日志和停止命令。主配置指定后台运行。重复启动不会派生第二个相同输出目录的调度器。可用`--once`取得单次快照。

患者步骤显示完成、失败、进行中、待完成、待核实、不适用及总人数；区域拟合和bootstrap显示实际任务数。原始计划名单不因失败缩小。心跳不等于算法进展，未知状态不当作失败或待完成。

监测输出在本run的monitor目录，含JSON、TSV及离线HTML。浏览器打开HTML不需要服务进程。截图应包括run_id、更新时间、阶段名称和错误摘要。

## 安全停止和恢复

`bash scripts/stop_project.sh --config config/project.yaml`只操作已核验的本项目调度器。默认停止新派发，当前任务继续完成；内部验证在R子任务检查点处可安全返回。超时只报告仍在收尾，不自动强杀。

再次启动将校验完成产物并续跑。若需要重试失败患者，把`runtime.yaml`中的`retry_failed_on_resume`设为true；成功患者行被保留，仅重试失败者，并使依赖这些结果的下游缓存失效。修复后可恢复false。

更改CPU预算或刷新间隔可兼容续跑；更改科学参数、输入、代码或实际软件版本会创建新run。不要手工删除锁或编辑finished.flag。状态损坏会要求核查而不是自动猜测。

## 截图交互与记录

每次反馈附run_id、监测截图和具体结果编号。使用`ANALYSIS_RECORD_TEMPLATE.md`记录观察、解释、决定及待验证事项。只有截图时，结论标为“依据截图”；拿到CSV/JSON复核后再升级验证状态。保留不显著、不可估计和未执行结果。

## Git与物理传递

代码、配置模板、4个公开图谱及验收摘要通过`https://github.com/zhenzonglin/WMH-GMV-PCA`传递。患者表、原始影像、拟合模型和完整结果放在配置指向的私有存储；物理包接收后校验SHA256。安装目录`.venv/.Rlib`不直接作为跨系统可移植环境；优先按锁文件重建。

完整记录备份解压到独立历史目录，仅供核查；旧run包含原绝对路径与主机身份，不能当作目标工作站正在运行的任务。可单独复制其私有配置作为起点，再审计实际数据。详见TRANSFER.md。
