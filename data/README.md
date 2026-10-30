# 私有输入目录

此目录仅README进入Git。默认名册路径是`data/LS/participants.tsv`、`data/C3/participants.tsv`和`data/GB/participants.tsv`；对应衍生数据目录为各队列下的`derivatives/`。这些是可修改的路径模板，不表示真实数据已经放入。

用config/paths.yaml或config/local中的私有配置指向实际名册、WMH20、GM119、真实ICV、总量、临床表和病灶掩膜。ID连接和字段编码见config/variables.yaml及docs/CLINICAL_CODEBOOK.md。

真实患者表、影像、病灶和含患者信息的派生结果通过机构批准的私有存储传递。当前本地工作区只有合成输入；真实三队列资料在远程工作站，未包含在本地物理包中。
