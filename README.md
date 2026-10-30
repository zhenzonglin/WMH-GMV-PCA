# WMH-GMV-PCA

WMH与GMV分别建立健康偏离、分别提取空间模式，再形成连续联合表型。LS为健康参考，C3为发现/建模队列，GB为冻结表征验证队列。正式年龄50–75岁。

本项目读取既有WMH20、原始GM119、真实ICV及急性病灶掩膜，不重做WMH校正，不剔除病灶累及GMV，不做中心校正，不纳入中心或MRI距发病时间。条件失配、脑龄和SuStaIn不在本项目中。

## 首次使用

1. 克隆`https://github.com/zhenzonglin/WMH-GMV-PCA.git`。4个公开图谱、标签和许可已随Git传递，在根目录执行`sha256sum -c resources/atlases/SHA256SUMS`校验。
2. 安装Linux/WSL环境、`uv`、R和基本编译工具；在`config/paths.yaml`配置实际路径。依赖以`requirements.lock`及`renv.lock`为准，不能只安装名称相同的软件包。
3. 执行`bash scripts/setup_project.sh --config config/project.yaml`安装项目隔离环境；R包写入配置的项目库，不修改输入数据。
4. 执行`.venv/bin/python scripts/audit_inputs.py --config config/project.yaml`先检查远程数据和临床变量。
5. 修改配置中的路径、字段编码和计算预算后启动。

```bash
bash scripts/start_project.sh --config config/project.yaml
.venv/bin/python scripts/monitor_project.py --config config/project.yaml
bash scripts/stop_project.sh --config config/project.yaml
```

单次查看使用`monitor_project.py --config config/project.yaml --once`。停止默认等待当前计算完成；脚本返回“draining”表示仍在收尾。再次手动调用启动命令会校验并继续兼容的运行。失败任务默认不自动重试。

## 配置

主入口`config/project.yaml`引用路径、分析、变量、运行四份配置。相对数据路径始终相对`project_root`。所有源ID按字符串读取。科学参数变化、源数据变化、代码或实际软件版本变化建立新运行，不覆盖旧结果。

`config/paths.yaml`使用项目相对路径模板，真实机器路径应放在私有配置中。原工作站路径保存在物理包的`config/local/workstation-original-20261030/`，仍需远程核查。临床表路径默认为空，不能由变量字典推断已有患者数据。JSON、TSV及已有轻量原始特征包均有读取适配。

`config/variables.yaml`中的性别、是/否、缺失码必须按实际源表核对。1/2编码不能未经映射直接作为0/1使用。无记录不会自动解释为无事件。

## 文档与输出

- `docs/statistical_analysis_plan.md`：完整统计方法、预期结果及解释边界。
- `docs/REMOTE_RUN.md`：远程安装、启动、监测、停止、恢复与截图交互。
- `docs/IMAGING_METHODS.md`、`docs/CLINICAL_METHODS.md`：实现和方法细节。
- `docs/ANALYSIS_RECORD_TEMPLATE.md`：每轮结果讨论与决定记录。
- `docs/VERIFICATION.md`：本地软件验收及未完成的远程验证。
- `docs/TRANSFER.md`：代码包、图谱包与SHA256校验步骤。

每个run保存输入清单、配置、软件和代码身份、任务日志、失败原因和产物哈希。监测分开显示患者步骤与区域/模型/bootstrap任务。PCA为0维是允许结果，不会为了获得表型而更改阈值。

每队列的`statistician_package`包含影像特征宽表、随访长表、事件和用药快照、变量说明、缺项清单与哈希。B3/B4不拟合出血模型。用药快照不代表连续用药区间。

## 合成软件验证

```bash
.venv/bin/python scripts/make_synthetic.py --config config/synthetic_v2.yaml
bash scripts/start_project.sh --config config/local/smoke_v2/project.yaml
.venv/bin/python -m pytest -q
```

合成输入、简化维数/重复次数以及验证输出均与正式配置分开。首次生成后若要重建，使用明确设置`overwrite: true`的生成配置，不清空真实数据目录。代码测试通过不代表远程真实队列分析已经完成。

`synthetic_v2.yaml`为正常可识别场景；`synthetic.yaml`保留早期弱识别场景，用于研究数值诊断，不能把它的历史失败改写为成功。零个稳定PC是允许结果，软件必须明确输出未增加维度而非拼出空交互项。

公开图谱及标签已随仓库提供。仅需重新准备资源时，配置FSL安装位置后使用`.venv/bin/python scripts/prepare_atlases.py --config config/project.yaml --sources config/atlas_sources.yaml`；详见`docs/ATLAS_RESOURCES.md`。

## 传递

Git传代码、配置模板、文档、依赖锁、4个审核后的公开图谱、原XML和合成验收摘要。患者数据、模型、日志、私有配置、本机环境及完整历史由`.gitignore`排除。公开资源使用明确文件名单，不对患者NIfTI开放通配上传。

工作站有网络并具备R/编译依赖时，仅克隆仓库即可重建软件和获取项目图谱；真实患者输入另行配置。本机记录和环境分别提供物理备份，均不是首次运行的必需包。分类、校验与迁移边界见`docs/TRANSFER.md`。
