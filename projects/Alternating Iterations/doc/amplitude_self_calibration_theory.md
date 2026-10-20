# 联合幅度自校准理论说明

## 1. 文档目的

本文说明 `Alternating Iterations` 项目中幅度自校准方法的数学模型、首次更新能够改善全带输出恢复的原因，以及高阶算子重复迭代后可能退化的原因。

当前实验假设未知反馈滤波器为线性相位 FIR。完成公共时延和常相位对齐后，主要剩余失配表现为通带幅度不平坦。

需要特别说明：当前实验采用 `R=250`，对应 501 抽头。因此该方法不应再称为“低维算子估计”，更准确的名称是“高阶结构化幅度算子估计”。

## 2. 反馈模型与传统覆盖模型

设全带 BLMP 基函数矩阵为

\[
\mathbf y_0=\mathbf\Psi(\mathbf x)\boldsymbol\theta_0,
\]

其中 `x` 为 PA 输入，`y_0` 为 PA 全带输出，\(\boldsymbol\theta_0\) 为 PA 参数。对齐后，真实反馈信号可以表示为

\[
\mathbf y_{\mathrm{bf}}
 =\mathbf H_2\mathbf A_\star
 \mathbf\Psi(\mathbf x)\boldsymbol\theta_0+\mathbf n,
\]

其中 \(\mathbf H_2\) 是已知数字覆盖滤波器，\(\mathbf A_\star\) 是未知滤波器 \(H_1\) 的相对幅度响应，\(\mathbf n\) 表示噪声和未建模误差。

传统覆盖方法忽略 \(\mathbf A_\star\)，使用

\[
\mathbf y_{\mathrm{bf}}
 \approx\mathbf H_2\mathbf\Psi(\mathbf x)\boldsymbol\theta.
\]

因此，未知滤波器的频率相关幅度会被错误地吸收到 PA 参数中。若记 \(\mathbf\Psi_2=\mathbf H_2\mathbf\Psi\)，则传统最小二乘估计的误差可近似写为

\[
\widehat{\boldsymbol\theta}^{(0)}-\boldsymbol\theta_0
\approx
(\mathbf\Psi_2^H\mathbf\Psi_2)^{-1}
\mathbf\Psi_2^H\mathbf H_2
(\mathbf A_\star-\mathbf I)
\mathbf\Psi\boldsymbol\theta_0.
\]

这说明传统方法的主要问题是回归矩阵失配，而不只是一个整体复增益误差。由于最终恢复使用未滤波的全带基函数，该失配会导致 PA 全带输出恢复偏差。

## 3. 当前幅度算子的实际形式

程序中的 `estimate_amplitude_operator.m` 使用

\[
A(q)=a_0+\sum_{r=1}^{R}a_r(q^r+q^{-r}).
\]

对应回归矩阵的第 `r+1` 列为

\[
z_r(n)=z_0(n-r)+z_0(n+r),
\]

其中 \(z_0=H_2\widehat y_0\)。因此，当前方法具有以下特征：

- 作用范围是双边 \([-R,R]\)；
- 由于正负抽头被强制相等，只估计 \(R+1\) 个独立实系数；
- 最终因果 FIR 为
  \[
  h_A=[a_R,\ldots,a_1,a_0,a_1,\ldots,a_R]^T;
  \]
- 因果实现引入 `R` 个样点群时延；
- 频率响应为
  \[
  A(e^{j\omega})=a_0+2\sum_{r=1}^{R}a_r\cos(r\omega),
  \]
  因而是实数、频率偶对称的零相位幅度响应。

在当前设置 `R=250` 时，算子长度为 501，独立参数数量为 251。它相对于一般单边复 FIR 仍有结构约束，但自由度已经较高，不属于原计划中 `R=3,5,8,10,20` 意义下的低维模型。

## 4. 单边估计与双边估计

单边 FIR 可写为

\[
A_{\mathrm{one}}(q)=\sum_{r=0}^{R}b_rq^{-r}.
\]

它只使用当前和过去样点，能够描述一般因果 FIR，但幅度和相位会同时变化，也不保证线性相位。

当前双边对称形式为

\[
A_{\mathrm{two}}(q)=a_0+\sum_{r=1}^{R}a_r(q^r+q^{-r}).
\]

它在完成时延和常相位对齐后只描述相对幅度响应，与“未知滤波器仍为线性相位”的假设一致。

因此，当前代码是“算子双边、参数单边”：完整算子为双边对称 FIR，但回归只估计一侧的独立系数。

如果使用一般单边 FIR，虽然拟合能力更强，但它可能把剩余时延、相位误差、PA 模型误差和对齐误差一起吸收到算子中，物理可解释性和可辨识性会下降。

## 5. 首次自校准为什么能够改善恢复

传统覆盖方法先得到粗略恢复

\[
\widehat y_0^{(0)}=\mathbf\Psi\widehat{\boldsymbol\theta}^{(0)}.
\]

随后通过

\[
\widehat{\mathbf A}^{(1)}
=\arg\min_{\mathbf A}
\left\|
\mathbf y_{\mathrm{bf}}
-\mathbf H_2\mathbf A\widehat y_0^{(0)}
\right\|_2^2
 +\lambda_a\|\mathbf a\|_2^2
\]

估计幅度算子。

若

\[
\widehat y_0^{(0)}=y_0+e_0,
\]

则在频域中有近似关系

\[
\widehat A^{(1)}(\omega)
\approx
H_1(\omega)
\frac{Y_0(\omega)}{\widehat Y_0^{(0)}(\omega)}.
\]

当粗略输出已经保留真实输出的主要频谱和非线性结构时，误差比值接近 1，于是

\[
\widehat A^{(1)}(\omega)\approx H_1(\omega).
\]

即使粗略输出并不完全准确，估计结果仍能形成一个有效的等效补偿算子。将其施加到全部 BLMP 基函数后，第二轮回归矩阵

\[
\mathbf H_2\widehat{\mathbf A}^{(1)}\mathbf\Psi
\]

比传统矩阵 \(\mathbf H_2\mathbf\Psi\) 更接近实际反馈链路对应的回归空间，从而降低滤波器失配向 PA 参数估计中的泄漏。

因此，自校准改善的本质是：

> 利用粗略 PA 输出恢复反馈链路中缺失的频率相关幅度结构，再用修正后的回归矩阵重新估计 PA 参数。

## 6. 为什么算子不一定等于真实 `H1`

由于

\[
(c\mathbf A)(\boldsymbol\theta/c)
=\mathbf A\boldsymbol\theta,
\]

幅度算子与 PA 参数之间存在尺度歧义。程序使用

\[
A(0)=a_0+2\sum_{r=1}^{R}a_r=1
\]

进行归一化，固定算子的中心频率增益。因此，估计结果主要表示未知滤波器的相对幅度形状，整体增益由 PA 参数承担。

此外，因算子输入是带误差的 \(\widehat y_0^{(0)}\)，估计结果还可能包含初始 PA 恢复误差的补偿成分。因此更准确的解释是

\[
\widehat A^{(1)}
=A_\star+\Delta A_{\mathrm{PA}}
+\Delta A_{\mathrm{noise}}.
\]

当前方法验证的是“等效幅度校正能否改善全带恢复”，而不是仅验证物理滤波器 `H1` 的精确辨识。

## 7. 为什么重复迭代可能越来越差

交替过程可写为

\[
\widehat{\boldsymbol\theta}^{(i)}
\rightarrow
\widehat{\mathbf A}^{(i+1)}
\rightarrow
\widehat{\boldsymbol\theta}^{(i+1)}.
\]

它近似求解联合目标

\[
\min_{\boldsymbol\theta,\mathbf a}
\left\|
\mathbf y_{\mathrm{bf}}
-\mathbf H_2\mathbf A(\mathbf a)
 \mathbf\Psi\boldsymbol\theta
\right\|_2^2.
\]

反馈拟合误差可能下降，但真实输出 NMSE 不一定下降，原因包括：

1. `R=250` 带来较高算子自由度，算子可以吸收 PA 模型误差；
2. 算子估计使用上一轮带误差的恢复输出，属于带误差回归；
3. 算子和 PA 参数存在乘积型可辨识性问题；
4. 反馈频带之外的信息无法由带限数据直接约束；
5. 全量更新会把上一轮误差直接带入下一轮。

因此，首次更新通常先修正主要幅度失配，后续更新则可能开始过拟合恢复误差。这解释了当前实验中“一次迭代最好”的现象。

## 8. 理论结论

当前方法应表述为：

> 一种基于实系数、偶对称和岭正则约束的高阶线性相位等效幅度算子交替自校准方法。

在粗略 PA 输出仍保留主要真实输出结构、反馈链路满足线性相位假设且信号在目标频带内具有充分激励时，首次算子更新可以降低模型侧回归矩阵与真实反馈链路之间的失配，从而改善 PA 全带输出恢复。由于算子阶数较高，该方法没有输出 NMSE 单调下降的理论保证，实际使用应结合早停、阻尼或更强的平滑正则化。
