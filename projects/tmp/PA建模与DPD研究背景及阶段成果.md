# PA 建模与 DPD 研究背景及阶段成果

> 阶段性概述，依据当前 MATLAB 工作区中的第三章 LD-PLS 代码、第四章幅度自校准理论文档及 BW=80 MHz 仿真结果整理。本文用于梳理研究脉络，不代替论文正文。

## 一、研究背景与总体问题

射频功率放大器（PA）在接近饱和区工作时具有非线性和记忆效应，会使发射信号产生带内失真与带外频谱扩展。数字预失真（DPD）需要较准确的 PA 行为模型，因此 PA 参数辨识和输出重建是后续线性化的基础。

实际系统中的反馈接收链路通常只提供有限带宽的观测。PA 输出的一部分频谱信息被滤除后，依据反馈数据构造的记忆多项式基函数可能高度相关，参数估计变得不稳定；若模型侧使用的滤波器与真实反馈链路不一致，还会产生额外的结构化模型失配。本研究围绕这两个相继出现的问题展开：第三章研究**已知且匹配的反馈滤波器条件下，如何稳定辨识**；第四章研究**未知带内幅度起伏导致滤波器不匹配时，如何自校准模型侧观测算子**。

设 $\mathbf\Psi$ 为由 PA 输入构造的全带记忆多项式基函数矩阵，$\boldsymbol\theta$ 为 PA 模型参数，$\mathbf\Phi$ 为经过反馈链路滤波后的带限基函数矩阵。研究的基本关系为

$$
\mathbf y\approx\mathbf\Psi\boldsymbol\theta,\qquad
\mathbf y_{\mathrm{bf}}\approx\mathbf\Phi\boldsymbol\theta,\qquad
\widehat{\mathbf y}=\mathbf\Psi\widehat{\boldsymbol\theta}.
$$

其中 $\mathbf y_{\mathrm{bf}}$ 是用于辨识的带限反馈观测，$\widehat{\mathbf y}$ 是重建的 PA 全带输出。最终重建始终使用未经过反馈滤波的 $\mathbf\Psi$。

## 二、第三章：带限反馈下的 LD-PLS 辨识

### 2.1 研究问题

在反馈带宽受限、模型侧滤波器与实际反馈链路匹配的前提下，带限记忆多项式（BLMP）能够把反馈滤波作用纳入基函数构造。但带限基函数之间的强相关会使直接估计 PA 参数变得不稳定。第三章的工作重点是改善这一辨识过程，并由估计出的参数重建全带 PA 输出。

### 2.2 已实现的方法

现有 LD-PLS 实现以记忆多项式基函数为基础，主要包含以下步骤：

1. 根据 PA 输入分别构造全带基函数 $\mathbf\Psi$ 和带限基函数 $\mathbf\Phi$；
2. 选取线性主导基函数，先从带限反馈观测中估计其系数；
3. 将反馈观测和其余基函数中与该线性分量相关的部分投影剔除；
4. 对非线性残差应用偏最小二乘（PLS）与正则化，估计剩余模型参数；
5. 使用全带基函数重建 PA 输出，并与 MG、MR、SBL 等方法作对照。

这一方法针对的是**滤波器匹配条件下的病态参数辨识**。仓库已包含 LD-PLS 函数、PLS 求解函数及对照实验脚本。现有第三章项目目录没有独立的最终数值结果表，因此本文不对第三章的 NMSE 改善幅度作定量声称。

### 2.3 第三章的理论依据

令 $\mathbf H$ 表示已知且与实际反馈链路匹配的带限滤波算子，则 $\mathbf\Phi=\mathbf H\mathbf\Psi$，带限观测可写为

$$
\mathbf y_{\mathrm{bf}}=\mathbf\Phi\boldsymbol\theta_\star+\boldsymbol\eta.
$$

如果直接采用最小二乘，且 $\mathbf\Phi$ 满列秩，则参数误差为 $\widehat{\boldsymbol\theta}_{\mathrm{LS}}-\boldsymbol\theta_\star=\mathbf\Phi^\dagger\boldsymbol\eta$，从而满足

$$
\|\widehat{\boldsymbol\theta}_{\mathrm{LS}}-\boldsymbol\theta_\star\|_2
\leq \frac{\|\boldsymbol\eta\|_2}{\sigma_{\min}(\mathbf\Phi)}.
$$

当带限处理使某些基函数列接近线性相关时，$\sigma_{\min}(\mathbf\Phi)$ 变小，噪声和模型误差可能被放大。这是第三章采用投影残差化与潜变量回归的动机；上式本身并不证明 LD-PLS 在所有数据上优于其他辨识器。

现有实现将 $\mathbf\Phi=[\boldsymbol\phi_{\mathrm L},\mathbf\Phi_{\mathrm R}]$，其中 $\boldsymbol\phi_{\mathrm L}$ 是选定的线性主导基函数，并先计算

$$
\widehat\theta_{\mathrm L}
=\frac{\boldsymbol\phi_{\mathrm L}^{H}\mathbf y_{\mathrm{bf}}}
{\boldsymbol\phi_{\mathrm L}^{H}\boldsymbol\phi_{\mathrm L}},
\qquad
\mathbf P_{\mathrm L}
=\boldsymbol\phi_{\mathrm L}
(\boldsymbol\phi_{\mathrm L}^{H}\boldsymbol\phi_{\mathrm L})^{-1}
\boldsymbol\phi_{\mathrm L}^{H}.
$$

随后构造残差观测 $\mathbf r=(\mathbf I-\mathbf P_{\mathrm L})\mathbf y_{\mathrm{bf}}$ 与残差基函数 $\mathbf Q=(\mathbf I-\mathbf P_{\mathrm L})\mathbf\Phi_{\mathrm R}$。PLS 从 $\mathbf Q$ 中抽取少量与 $\mathbf r$ 相关的潜变量，再在潜变量空间使用岭正则估计非线性部分。对标准化后的潜变量矩阵 $\mathbf T$，其回归步骤可概括为

$$
\widehat{\mathbf c}
=(\mathbf T^H\mathbf T+\lambda_\theta\mathbf I)^{-1}
\mathbf T^H\mathbf r.
$$

这里的 $\mathbf T$ 和 $\mathbf r$ 应理解为与代码相应的标准化量；实际系数还需经过潜变量变换与尺度还原。现有代码按预设比例 $S$ 取潜变量个数，不能仅凭理论稿中的表述声称已实现自适应选择。

## 三、当前成果：非平坦反馈链路的结构化幅度自校准

### 3.1 从第三章延伸出的新问题

第三章要求模型侧带限算子与真实反馈链路相匹配。当模拟反馈滤波器存在未知的带内幅度起伏时，仅用已知数字滤波器构造 $\mathbf\Phi$，会使模型侧基函数矩阵与实际观测不一致。即使参数辨识器能够处理基函数相关性，也不能自动消除这种算子失配。

当前程序针对模拟的反馈链路幅度失配。对**待估计的残余 FIR 算子，只施加实系数、偶对称这一结构**；程序并未施加频响非负约束，也未估计一般复数频响。算子的增益归一化、信号对齐和正则化均按现有代码执行。

### 3.2 第一个滤波器 H1 的实验设置

主实验用两个串联 FIR 模拟反馈链路：PA 全带输出 $y_0$ 先经过第一个滤波器 **H1**，再经过已知数字滤波器 **H2**。H1 用于产生待校准的带内起伏；H2 是模型侧传统覆盖时已知的滤波器。H1 的系数虽然在仿真程序中可得，但在 Initial 和自校准组中**不作为已知覆盖滤波器使用**；只有 Ref 组使用真实的 $h_1*h_2$。

为使扫描时 H1 的低通部分保持一致，程序先在循环外设计一次平坦低通 $h_{\mathrm{flat}}$：采样率 $F_s=1000$ MHz，$BW=80$ MHz，通带边缘 $BW+5=85$ MHz，阻带起点 $BW+10=90$ MHz；设计通带纹波为 0.01 dB，阻带幅度误差为 $10^{-4}$。阶数由 `firpmord` 确定，系数由 `firpm` 生成，随后各纹波点共用该 $h_{\mathrm{flat}}$。

每个纹波点用一个对称三非零抽头滤波器叠加余弦起伏。设标称峰峰纹波为 $R_p$ dB，程序采用

$$
\varepsilon=\frac{10^{R_p/20}-1}{10^{R_p/20}+1},
\qquad
\tau=\operatorname{round}\!\left(\frac{n_{\mathrm{Ripple}}F_s}{BW+5}\right),
\qquad n_{\mathrm{Ripple}}=5.
$$

当前参数下 $\tau=59$ 个样点。$h_{\mathrm{ripple}}$ 在第 $0,\tau,2\tau$ 个抽头处分别取 $\varepsilon/2,1,\varepsilon/2$，其余抽头为零。H1 由

$$
h_1=h_{\mathrm{flat}}*h_{\mathrm{ripple}}
$$

构成；忽略纯时延时，纹波部分的频率响应为 $1+\varepsilon\cos(2\pi f\tau/F_s)$。这样的设置使低通部分在整个扫描中固定，而余弦纹波幅度随 $R_p$ 调整。标称纹波由上述公式设置，**实际纹波**则按程序用 `freqz` 在 $0\leq f\leq BW+5$ MHz 内计算 $20\log_{10}(\max|H_1|/\min|H_1|)$；因此二者会有小幅差别，例如标称 1.9 dB 对应实测约 1.909 dB。

H2 在程序中使用 $BW=80$ MHz 的通带边缘、$BW+3=83$ MHz 的阻带起点，通带幅度误差 `Dpass2 = 0.0057501127785`、阻带幅度误差 `Dstop2 = 1e-4`，由 `firpmord` 和 `firpm` 设计。在每个扫描点，观测按 $y_f=H1(y_0)$、信号对齐、$y_{\mathrm{bf}}=H2(y_f)$ 的顺序生成。

### 3.3 算法流程与理论成果

幅度自校准先沿用第三章的传统覆盖，获得粗略 PA 全带输出；再结合带限反馈观测，通过带正则项的回归估计残余幅度算子；随后将该算子并入模型侧滤波器，重新构造带限基函数，并再次使用 LD-PLS 估计 PA 参数。当前主实验采用**一次自校准更新**。

已有理论文档分析了传统覆盖的回归矩阵失配和算子修正后的参数误差来源。以下推导用于解释程序结果；具体估计步骤以 `estimate_amplitude_operator.m` 和 `run_amplitude_alternating_iterations.m` 为准。反复交替更新的反馈拟合改善不必然转化为全带 NMSE 的单调改善。

### 3.4 当前自校准方法的理论依据

令 $\mathbf H_{\mathrm L}$ 表示已知数字覆盖滤波器 H2，$\mathbf A_\star$ 表示对齐后 H1 对观测造成的待校准部分。在忽略或并入模型误差的写法下，反馈观测为

$$
\mathbf y_{\mathrm{bf}}
=\mathbf H_{\mathrm L}\mathbf A_\star\mathbf\Psi\boldsymbol\theta_\star
+\boldsymbol\eta.
$$

传统覆盖令 $\mathbf A_\star\approx\mathbf I$，对应 $\mathbf\Phi_0=\mathbf H_{\mathrm L}\mathbf\Psi$。若 $\mathbf\Phi_0$ 满列秩，并以伪逆估计参数，则误差可以精确分解为

$$
\widehat{\boldsymbol\theta}^{(0)}-\boldsymbol\theta_\star
=\mathbf\Phi_0^\dagger\mathbf H_{\mathrm L}
(\mathbf A_\star-\mathbf I)\mathbf\Psi\boldsymbol\theta_\star
+\mathbf\Phi_0^\dagger\boldsymbol\eta.
$$

第一项是未知幅度响应造成的**结构化算子失配**，第二项是噪声及并入其中的模型误差。LD-PLS 能缓解病态辨识，但仅更换求解器不能使第一项自动消失。这也是第四章需要先修正模型侧覆盖算子的原因。

当前程序估计一个实系数、偶对称的残余算子，写作

$$
\mathbf A(\mathbf a)
=a_0\mathbf I+\sum_{r=1}^{R}a_r(\mathbf T_r+\mathbf T_{-r}),
\qquad
A(e^{j\omega})=a_0+2\sum_{r=1}^{R}a_r\cos(r\omega),
$$

其中 $\mathbf T_r$ 为时移算子。实偶对称结构使居中表示的频率响应为实偶函数；程序生成因果对称 FIR 时相应引入固定时延。实验取 $R=250$，对应长度为 501 的 FIR、251 个独立实参数，因此它是**高阶结构化算子**。实偶对称本身并不保证频响处处为正。

给定初始全带输出 $\widehat{\mathbf y}^{(0)}$，令 $\mathbf Z(\widehat{\mathbf y}^{(0)})$ 的第零列为 $\mathbf H_{\mathrm L}\widehat{\mathbf y}^{(0)}$，其余列为该信号相应的正、负对称移位之和，则

$$
\mathbf H_{\mathrm L}\mathbf A(\mathbf a)\widehat{\mathbf y}^{(0)}
=\mathbf Z(\widehat{\mathbf y}^{(0)})\mathbf a.
$$

程序先对粗略全带输出应用 H2，再用 `circshift` 构造 $\mathbf Z$ 的各对称移位列；只对反馈目标进行一次公共时延与常相位对齐，并裁掉两端的滤波暂态和环绕样本。设裁剪后的矩阵和目标为 $\mathbf Z_v,\mathbf y_v$，实际求解的是实参数岭回归

$$
\mathbf a_{\mathrm{raw}}
=\bigl[\operatorname{Re}(\mathbf Z_v^H\mathbf Z_v)+\lambda_a\mathbf I\bigr]^{-1}
\operatorname{Re}(\mathbf Z_v^H\mathbf y_v),
\qquad
\lambda_a=\eta\,\frac{\operatorname{tr}(\operatorname{Re}(\mathbf Z_v^H\mathbf Z_v))}{R+1},
$$

其中 `operator_eta = 1e-4`。随后按 $g=a_{\mathrm{raw},0}+2\sum_{r=1}^{R}a_{\mathrm{raw},r}$ 计算中心频率增益，并令 $\mathbf a=\mathbf a_{\mathrm{raw}}/g$；再把 $\mathbf a$ 排列为因果对称 FIR $h_A$，构造 $h_{\mathrm{model}}=h_A*h_2$。这个**先岭回归、后归一化**的步骤用于处理算子与 PA 参数之间的尺度歧义；它不同于理论稿中直接施加等式约束的回归。

首次校准的精度仍取决于初始输出是否足够准确。若 $\widehat{\mathbf y}^{(0)}=\mathbf y_{\mathrm m}+\mathbf e_0$，在理想线性观测模型下，真实算子处的回归残差包含 $\boldsymbol\eta-\mathbf H_{\mathrm L}\mathbf A_\star\mathbf e_0$。初始重建误差因此会进入算子估计；这里仅解释误差来源，不将理论稿针对另一种约束估计器的扰动界套用于当前代码。

使用估计算子后，模型侧带限基函数变为 $\mathbf\Phi(\mathbf a)=\mathbf H_{\mathrm L}\mathbf A(\mathbf a)\mathbf\Psi$。在满列秩、伪逆估计的分析模型下，其参数误差由剩余算子失配 $\mathbf A_\star-\mathbf A(\mathbf a)$ 与噪声共同决定。减小失配可以降低相应的误差**上界**，前提是更新后的 $\mathbf\Phi(\mathbf a)$ 没有变得更病态；这不是对实际 LD-PLS 输出 NMSE 必然改善的证明。

带限观测还存在根本的信息边界：若某个非零参数扰动 $\Delta\boldsymbol\theta$ 满足 $\mathbf\Phi(\mathbf a)\Delta\boldsymbol\theta=\mathbf0$，但 $\mathbf\Psi\Delta\boldsymbol\theta\neq\mathbf0$，则两组参数在反馈域无法区分，却给出不同的全带输出。因而，即使交替更新降低反馈拟合误差，也不能仅据此推断全带重建误差单调下降。

### 3.5 已完成的仿真结果

主实验在 BW=80 MHz、采样率 1000 MHz 的设置下，扫描反馈模拟滤波器的标称通带纹波 0.1–1.9 dB。三组对照分别为：使用真实复合滤波器的参考组（Ref）、仅使用已知数字滤波器的未校准组（Initial），以及一次自校准后的方法组（Iteration）。

| 标称通带纹波 | Ref NMSE | Initial NMSE | 一次自校准 NMSE |
|---|---:|---:|---:|
| 0.1 dB | −38.47 dB | −38.32 dB | −38.35 dB |
| 1.9 dB | −38.45 dB | −32.45 dB | −36.18 dB |

在 1.9 dB 纹波处，一次自校准较未校准改善 **3.73 dB**，与 Ref 的差距由 **6.00 dB** 缩小至 **2.27 dB**。完整的十个扫描点见 `projects/Alternating Iterations/results/result.md`。这一结果支持：在当前仿真条件下，结构化幅度自校准能够减轻反馈链路带内幅度失配造成的全带重建性能下降。

## 四、两阶段成果的关系与当前边界

第三章解决“**算子已匹配，但带限基函数难以稳定辨识**”；当前成果解决“**算子本身因未知幅度响应而失配**”。自校准修正模型侧带限基函数后，LD-PLS 继续承担稳定估计 PA 参数的任务。两部分共同形成从带限反馈观测到 PA 全带输出重建的研究主线。

目前仍需严格区分已验证结果与开放问题：

- 已完成的是一次幅度自校准的理论分析与仿真验证；实测数据尚未接入该框架。
- 当前实偶算子没有频域非负约束，不能保证其在所有情况下都可解释为严格的物理幅度响应。
- H1 的当前仿真实现由对称 FIR 构成，具有固定群延迟；程序没有估计一般的未知非线性相位或反馈带外响应。
- 现有结果不能推出任意迭代次数下全带 NMSE 单调下降，也不能直接视为 DPD 闭环线性化效果。

下一阶段拟在上述研究主线上寻找具有明确理论命题、可通过实验检验的新问题；具体方向尚未确定。

## 资料依据

- `projects/LD-PLS/main.m`、`projects/LD-PLS/PLS_Algorithms/get_OPLS.m`、`pls_B.m`：第三章实现与对照入口。
- `projects/Alternating Iterations/README.md`：当前方法、参数和实验说明。
- `projects/Alternating Iterations/doc/幅度自校准理论_统一符号与严格论证 (2).md`：理论分析及结论边界。
- `projects/Alternating Iterations/results/result.md`：BW=80 MHz 纹波扫描数据。
