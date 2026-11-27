# 仓库协作规范

## 项目结构与模块组织

本仓库是用于射频功率放大器建模和数字预失真（DPD）实验的 MATLAB 工作区。

- `common/` 存放可复用代码。功率放大器模型放在 `models/`，信号处理和指标计算工具放在 `tools/`，通用实验工具放在 `experiments/`。
- `projects/` 存放各算法项目，包括 `LD-PLS`、`MG`、`MR`、`SBL` 和 `MP_DPD`。项目入口脚本与项目专用辅助代码放在同一项目下。
- `datasets/` 存放共享输入数据。项目专用的小型 `.mat` 文件可以放在该项目的 `data/` 目录。
- `.matlab/startup.m` 负责初始化搜索路径。仓库根目录的 `startup.m` 是 MATLAB 可发现的入口，应保持轻量。

在 VS Code 中打开仓库根目录，使 MATLAB 从正确的目录启动，并将工作区及其子目录加入搜索路径。

## 构建、测试与开发命令

本仓库无需编译。请从仓库根目录执行以下命令。

```powershell
matlab -batch "startup; which get_NMSE"
matlab -batch "checkcode('common/tools/metrics/get_NMSE.m')"
matlab -batch "results=runtests('tests'); assertSuccess(results)"
```

第一条命令用于验证启动和搜索路径。第二条命令用 `checkcode` 对 MATLAB 代码做静态分析。存在 `tests/` 目录时，第三条命令可以运行测试套件。使用 `run('projects/<name>/main.m')` 运行实验前，先确认所需的数据和工具箱已就绪。

## 代码风格与命名规范

使用四个空格缩进，每行一条语句，二元运算符两侧留空格。函数名须与文件名一致。沿用 `get_`、`plot_` 和 `create_` 等描述性前缀。`main.m` 仅用作项目入口脚本。使用 `%` 注释说明算法假设和单位，避免逐行复述代码。路径通过 `fullfile` 相对项目或脚本的位置拼接，避免新增硬编码的本机路径。代码以简洁为主，仅在存在明确风险或公共接口确有需要时增加输入校验和异常分支。

MATLAB 绘图时，要用 `'Color', 'w'` 显式将图窗和坐标区背景都设为白色。对每个图窗和坐标区分别设置颜色，使白底显示不受桌面主题影响。

## 测试规范

MATLAB 单元测试放在对应项目的 `tests/` 目录，例如 `projects/Alternating Iterations/tests/`。测试文件使用 `test_<subject>.m` 命名，测试方法按验证的行为命名，例如 `test_identicalSignalsHaveZeroNMSE`。使用确定的输入和固定的随机种子，覆盖正常行为、数据维度、复数信号和无效参数。生成的结果放在同一项目的 `results/` 目录，测试期间须保留已提交的 `.mat` 测试数据原样。

测试脚本统一放在各自项目的 `tests/` 目录，其中包括 `test_*.m` 单元测试和 `diag_*.m`、`verify_*.m` 诊断及验证脚本。

## 提交与拉取请求规范

提交标题和正文统一使用中文。代码标识符、文件路径和专有缩写可以保留原文。标题应简洁明确，例如 `新增：LD-PLS 正则化参数扫描` 或 `修复：参考功率为零时的 NMSE 计算`。每次提交保持聚焦。

拉取请求应说明算法改动、MATLAB 版本、所需工具箱、数据来源、验证命令和改动前后的数值结果。信号表现发生变化时，在拉取请求中附图。生成的图文件和大型数据集仅在明确要求时纳入提交。
