# matlab_Prime_Power_Tilings_tile_partitions
MATLAB codes for the paper "Tilings of Integer Sets of Prime-Power Cardinality"

# 单素数幂密铺枚举器（26-11-17 版本）

本程序依据论文 *Fundamental periods and complements of prime-power integer tiles*
（定理 3.1、4.2、推论 4.4、命题 4.9）与
*Prime_Power_Tilings_tile_partitions*（26-11-16），枚举基数 $|A|=p^\alpha$
（$p$ 素数）的整数密铺集 $A$ 在各可行模数 $\mathbb Z_M$ 下的**全部铺法**，
并按三种等价关系分类给出类代表。

---

## 1. 数学背景速览

设 $|A| = p^\alpha$，$A(x)=\sum_{a\in A}x^a$。

- **素幂层与纤维公约数**：$\Phi_{p^j}(x)\mid A(x)$ 当且仅当 $j\in\{\nu_1<\dots<\nu_k\}$；
  $\Delta_i=\gcd(\{a\bmod p^{\nu_i}\})$，记 $P=p^{\nu_k}$、$D=\Delta_k$。
- **基本周期**（定理 3.1）：恰为 $\{P\cdot d:\ d\mid D\}$；$P\mid M$ 是 $\mathbb Z_M$ 存在铺法的必要条件。
- **截面塔枚举**（定理 4.2）：以 $Pd$ 为周期的补集个数
  $F_A(Pd)=p^{\sum_i p^{\nu_i-i}\gcd(d,\Delta_i)}$，所有补集由逐层截面函数生成。
- **类数**（推论 4.4）：$E_A(Pd)=\sum_{e\mid d}\mu(d/e)F_A(Pe)$ 为基本周期恰为 $Pd$
  的补集个数；循环平移等价类数 $=E_A(Pd)/(Pd)$。
- **三种等价关系**：
  1. **循环平移等价**：$B\sim B+t\ (t\in\mathbb Z_M)$；
  2. **仿射等价**：$B\sim uB+t\ (u\in\mathbb Z_M^\times,\ t\in\mathbb Z_M)$。
     由 Tijdeman 膨胀引理，$\gcd(u,|B|)=1$ 时 $uB$ 仍是同一 $A$ 的补集；
  3. **砖块分拆等价**（equivalent as tile partitions）：
     $\mathcal P_M(B)=\{A+b : b\in B\}$ 作为无序砖块集合相同。
     **注意**：整体平移会移动全部砖块，故分拆等价**不**把平移视为等价——
     必须对全部铺法（每个平移类的 $M$ 个平移）计算分拆键。

---

## 2. 文件清单

| 文件 | 类型 | 说明 |
|---|---|---|
| `tiling_gui_v2.m` | 入口（GUI） | 图形界面。输入 A 与模数 M（留空 = 全部基周期），点击"枚举"。无绘图选项框。 |
| `tiling_report_v2.m` | 入口（命令行） | 主管线：`[txt, S] = tiling_report_v2(A, M, opts)`。生成完整报告文本与结构化结果，GUI 与演示脚本均调用它。 |
| `tile_enum_tilings.m` | 核心 | 枚举引擎（版本戳 `TILE_VER = '2026-11-17a'`）。计算算术数据、逐层截面塔枚举、平移类代表。无绘图、无"实际加载"调试行。 |
| `classify_tilings_v2.m` | 核心 | 仿射等价 + 砖块分拆等价分类：`out = classify_tilings_v2(A, M, reps, opts)`，输入为平移类代表集合。 |
| `reduce_tiling_period.m` | 辅助 | 约化到基周期：`(A, M) → (A*, M*)`。两种机制：A 为子群陪集之并（$F_A(\mathbb Z_M)=m^{\|B\|}F_A(\mathbb Z_d)$）；A 横截且 A mod d 可密铺 $\mathbb Z_d$。 |
| `tile_enum_tilings_test.m` | 测试 | 回归测试：论文 7 个算例 + 暴力枚举对拍 + 约化校验。全部通过。 |
| `demo_tiling_v2.m` | 演示 | 命令行演示三个算例（相当于 GUI 中输入后点"枚举"）。 |

依赖关系：`tiling_gui_v2 / demo_tiling_v2 → tiling_report_v2 → tile_enum_tilings + classify_tilings_v2 + reduce_tiling_period`。
将整个文件夹放入 MATLAB 路径即可运行。

---

## 3. 快速上手

```matlab
% GUI
>> tiling_gui_v2

% 命令行：全部基周期
>> txt = tiling_report_v2([0 3 12 15]);      % 基周期 {8, 24}
>> disp(txt)

% 命令行：只枚举 Z_24
>> txt = tiling_report_v2([0 3 12 15], 24);

% 命令行：约化到基周期
>> txt = tiling_report_v2([0 3 12 15], 24, struct('ReducePeriod', true));

% 回归测试
>> tile_enum_tilings_test
>> demo_tiling_v2
```

**选项 opts**（`tiling_report_v2` 与 GUI 共用）：

| 选项 | 默认 | 说明 |
|---|---|---|
| `ReducePeriod` | `false` | 勾选后先判定 A（或补集）是否为子群陪集之并，是则把 A 换成 A*、模数缩小后重新求解 |
| `MaxEnum` | `2e5` | 单层截面塔枚举上限，超出则只计数不列代表 |
| `MaxOps` | `5e7` | 仿射/分拆分类的运算量上限，超出则跳过该分类并告警 |
| `ShowReps` | `40` | 每个列表最多显示的类代表数 |

---

## 4. 输出格式

报告对每个模数 $M$ 依次给出（以 $A=[0\ 3\ 12\ 15]$、$M=24$ 为例）：

```
======== 模数 M = P*d = 24 ========
铺法总数 F_A(24) = 96，平移类 8，仿射类 5，分拆类 8
-- 平移等价类代表（8 类）--
    1)  B = [0 2 8 10 16 18]   （基本周期 8）
    2)  B = [0 1 2 6 7 8]      （基本周期 24）
    ...
-- 仿射等价类代表（5 类）--
    ...
-- 砖块分拆等价类（equivalent as tile partitions，共 8 类）--
    ...
```

**两条显示约定**（2026-11-17 修订）：

1. 所有输出的铺法 B 均已平移规范化使 $0\in B$（$B\mapsto B-\min(B)$；
   仍是同一平移类的合法铺法，不影响任何类数与计数）。
2. 同一基本周期层内的代表按**字典序**排列（`sortrows` 逐行比较）；
   层与层之间按基本周期从小到大排列。仿射、分拆类代表列表亦按字典序排列。

模数 M 可以比 A 的元素小：程序自动取 $A\bmod M$（横截时）求铺法，
如 $A=[0\ 3\ 10\ 13]$ 的基周期 $M=4<13$ 照常枚举。

---

## 5. 主要输出字段（`tile_enum_tilings`）

| 字段 | 含义 |
|---|---|
| `isTile` | A 是否为整数密铺集 |
| `p, k, nu, Delta` | 素数 p；素幂层数 k；素幂层 nu；纤维公约数 Delta |
| `P, D, PD` | 万能周期参数（基周期集 $=\{Pd: d\mid D\}$） |
| `table` | 每层 [d, Pd, F_A(Pd), E_A(Pd), 平移类数] |
| `numTilings / numClasses` | 铺法总数 / 平移类总数 |
| `representatives` | 平移类代表（cell，含 0、层内字典序） |
| `repsByLayer{t}` | 第 t 层（基本周期 $P\cdot\text{table}(t,1)$）的代表 |
| `report` | 报告文本 |

`classify_tilings_v2` 输出：

| 字段 | 含义 |
|---|---|
| `affine.count / reps / idx / canon` | 仿射类数 / 各类代表 / 每个平移代表所属类 / 各类典范形 $\min_{u,t}(uB+t)$ |
| `partition.count / reps / sizes` | 分拆类数 / 各类代表 / 各类实际铺法数（膨胀铺法的重复平移已去重） |
| `partitionSkipped` | 铺法总数超过 `MaxOps` 时为 `true`（跳过分拆分类） |

---

## 6. 验证情况

- **回归测试**（`tile_enum_tilings_test.m`）：论文 7 个算例（各层计数与类数逐一核对，
  含 $A=[0\ 10\ 20]$ 的 1978 平移类大算例）+ 13 组暴力精确覆盖对拍 + 约化校验，全部通过。
- **分拆类数公式**核对：$\#\{\text{分拆类}\}=p^{\sum_{i\le r}p^{\nu_i-i}d_i}$（论文记号），
  在 6 个算例上与程序逐一吻合，如 $A=[0\ 3\ 10\ 13]$ 的 $\mathbb Z_{20}$：64 铺法 → 2 分拆类
  （代表 $[0\ 2\ 4\ 6\ 8]$ 与 $[1\ 3\ 5\ 7\ 9]$，已用独立 Python 精确覆盖穷举仲裁）；
  $A=[0\ 3\ 12\ 15]$ 的 $\mathbb Z_{24}$：8 分拆类；$A=[0\ 10\ 20]$ 的 $\mathbb Z_{30}$：1 分拆类。
- **大算例性能**：$A=[0\ 10\ 20]$、$M=30$（59049 铺法）在 Octave 11.3 下约 51 秒
  （1978 平移类 / 300 仿射类 / 1 分拆类）。

兼容性：核心与报告在 MATLAB R2021a+ 与 GNU Octave 11 下均实测通过
（GUI 用 `uifigure/uigridlayout`，仅 MATLAB 可用）。
