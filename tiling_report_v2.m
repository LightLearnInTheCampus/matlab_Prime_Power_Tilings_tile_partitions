function [txt, S] = tiling_report_v2(A, Muser, opts)
%TILING_REPORT_V2  枚举密铺集 A 的全部基周期并给出三种等价类的代表
%
%   [txt, S] = tiling_report_v2(A)              % 枚举全部基周期
%   [txt, S] = tiling_report_v2(A, M)           % 只枚举 Z_M（M 须是 P 的倍数）
%   [txt, S] = tiling_report_v2(A, M, opts)
%
% 输出结构
%   1) 密铺集 A 的算术数据（p, k, nu, Delta, P, D）与全部基周期 {P*d : d | D}；
%   2) 对每个基周期 M：
%        - Z_M 中铺法总数、循环平移等价类数、仿射等价类数、砖块分拆等价类数；
%        - 平移等价类代表（每类一个，标注该铺法的基本周期）；
%        - 仿射等价类代表（B ~ u*B+t，gcd(u,M)=1）；
%        - 砖块分拆（equivalent as tile partitions）等价类代表。
%   A 中的元素可以比 M 大：程序自动用 A mod M 求铺法。
%
% 选项 opts
%   .ReducePeriod  约化到基周期（A 或补集为子群陪集之并时把 A 换成 A*、
%                  模数相应缩小）                       (默认 false)
%   .MaxEnum       单层枚举上限                          (默认 2e5)
%   .MaxOps        仿射分类运算量上限                    (默认 5e7)
%   .ShowReps      每个列表最多显示的类代表数            (默认 40)
%
% S 返回结构化结果（每个模数一个 S.periods{i} 条目）。

if nargin < 2 || isempty(Muser), Muser = []; end
if nargin < 3 || isempty(opts), opts = struct(); end
dflt = struct('ReducePeriod', false, 'MaxEnum', 2e5, 'MaxOps', 5e7, ...
              'ShowReps', 40);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = dflt.(fn{i}); end
end

L = {};
S = struct('input', A, 'reduced', false, 'periods', {{}});
cOpt = struct('MaxEnum', opts.MaxEnum, 'Report', false);

% ------------------------------------------------ 0. 第一次求解：算术数据
out0 = tile_enum_tilings(A, [], struct('MaxEnum', opts.MaxEnum, ...
             'Report', false, 'ReducePeriod', false));
if out0.n == 1
    % |A| = 1 平凡情形
    A0 = sort(double(A(:))' - min(double(A(:))'));
    txt = sprintf(['A = %s   （|A| = 1，平凡情形）\n' ...
        '唯一基周期 M = 1：唯一铺法 B = [0]（A mod 1 为单点），' ...
        '1 个平移类、1 个仿射类、1 个砖块分拆类。'], fmt_set(A0));
    S.periods{1} = struct('M', 1, 'trivial', true);
    return
end
if ~out0.isTile
    nn = out0.n;
    if nn > 1
        pp = factor(nn);  pp = pp(1);  kk = numel(factor(nn));
    else
        pp = 1;  kk = 1;
    end
    L{end+1} = '================ 密铺集 A 的算术数据 ================';
    L{end+1} = sprintf('A = %s   （|A| = %d = %d^%d）', fmt_set(out0.input), ...
        nn, pp, kk);
    L{end+1} = '结论：A 不是整数密铺集，无任何铺法（原因见下）。';
    L{end+1} = '';
    L{end+1} = out0.report;
    txt = strjoin(L, newline);
    return
end

% ------------------------------------------------ 0b. 可选：约化
Awork = out0.A0;                        % 规范化 A0（min = 0）
if opts.ReducePeriod
    if isempty(Muser), Mtarget = out0.PD; else, Mtarget = Muser; end
    [A2, M2, rinfo] = reduce_tiling_period(Awork, Mtarget, ...
                          struct('MaxEnum', opts.MaxEnum));
    if rinfo.reduced
        S.reduced = true;
        S.reduceInfo = rinfo;
        L{end+1} = '-------- 基周期约化：A（或其补集）是子群陪集之并 --------';
        for i = 1:numel(rinfo.steps)
            L{end+1} = sprintf('  第 %d 步：%s', i, rinfo.steps{i}); %#ok<AGROW>
        end
        L{end+1} = sprintf('  以下按约化后的 A* = %s（|A*| = %d）枚举。', ...
            fmt_set(A2), numel(A2));
        L{end+1} = '-----------------------------------------------------';
        Awork = A2;
        if ~isempty(Muser), Muser = M2; end
    end
end

% ------------------------------------------------ 1. Awork 的算术数据
outA = tile_enum_tilings(Awork, [], cOpt);
p = outA.p;  k = outA.k;  nu = outA.nu;  Delta = outA.Delta;
P = outA.P;  D = outA.D;

S.p = p; S.k = k; S.nu = nu; S.Delta = Delta; S.P = P; S.D = D;
S.Awork = Awork;

if isempty(Muser)
    dlist = all_divisors(D);
    basePeriods = P * dlist;
    S.basePeriods = basePeriods;
else
    basePeriods = Muser;
    dlist = NaN;
end

L{end+1} = '================ 密铺集 A 的算术数据 ================';
if S.reduced
    L{end+1} = sprintf('A = %s   （原始输入；约化后 A* = %s）', ...
        fmt_set(out0.input), fmt_set(Awork));
else
    L{end+1} = sprintf('A = %s   （规范化 A0 = %s，min = 0）', ...
        fmt_set(out0.input), fmt_set(Awork));
end
L{end+1} = sprintf('|A| = %d = %d^%d        p = %d,  alpha = k = %d', ...
    outA.n, p, k, p, k);
L{end+1} = sprintf('素幂层  nu_i = %s        P = p^{nu_k} = %d', ...
    fmt_set(nu), P);
L{end+1} = sprintf('fiber gcd  Delta_i = %s   D = Delta_k = %d', ...
    fmt_set(Delta), D);
L{end+1} = sprintf('典范等价密铺集  A* = %s', fmt_set(outA.canonicalTile));
if isempty(Muser) || isnan(dlist)
    L{end+1} = sprintf('全部基周期（Theorem 3.1）：M = P*d，d | D，即 M ∈ {%s}。', ...
        numjoin(basePeriods));
else
    L{end+1} = sprintf('指定模数 M = %d（须为 P = %d 的倍数）。', Muser, P);
end
L{end+1} = sprintf(['（求 Z_M 中的铺法时程序自动取 A mod M；' ...
    'M 比 A 的元素小时亦如此。）']);
L{end+1} = '';

S.periods = cell(1, numel(basePeriods));

% ------------------------------------------------ 2. 逐个模数枚举与分类
for t = 1:numel(basePeriods)
    Mp = basePeriods(t);
    out = tile_enum_tilings(Awork, Mp, cOpt);

    L{end+1} = repmat('=', 1, 64);
    if isempty(Muser)
        L{end+1} = sprintf('======== 基周期 M = P*d = %d（d = %d） ========', ...
            Mp, dlist(t));
    else
        L{end+1} = sprintf('======== 模数 M = %d ========', Mp);
    end

    if ~out.isTile || out.numTilings == 0
        L{end+1} = sprintf('Z_%d 中没有铺法（%s）。', Mp, ...
            strtrim(out.report));
        L{end+1} = '';
        S.periods{t} = struct('M', Mp, 'empty', true);
        continue
    end

    cls = classify_tilings_v2(Awork, Mp, out.representatives, ...
              struct('MaxOps', opts.MaxOps));

    L{end+1} = sprintf('A mod %d = %s（%d 个元素，无重合）', Mp, ...
        fmt_set(cls.Amod), cls.n);

    % ---- 计数总表
    L{end+1} = sprintf(['Z_%d 中铺法总数 F_A(%d) = %d；' ...
        '平移等价类 = %d；仿射等价类 = %s；砖块分拆等价类 = %d。'], ...
        Mp, Mp, out.numTilings, out.numClasses, ...
        class_str(cls.affine.count), cls.partition.count);

    % ---- 各层（按精确基本周期）信息
    L{end+1} = '  各层（按铺法的精确基本周期）：';
    L{end+1} = '    d      Pd     F_A(Pd)   E_A(Pd)  平移类数';
    for q = 1:size(out.table, 1)
        L{end+1} = sprintf('  %2d  %6d  %9d  %9d  %8d', out.table(q, :)); %#ok<AGROW>
    end

    % ---- 平移等价类代表
    L{end+1} = '';
    L{end+1} = sprintf(['-- 平移等价类代表（共 %d 类，模 %d；' ...
        '括号内为该铺法的基本周期）--'], out.numClasses, Mp);
    L = emit_reps(L, out.representatives, opts.ShowReps, layerPeriods(out));

    % ---- 仿射等价类
    L{end+1} = '';
    if strcmp(cls.skipped, 'affine')
        L{end+1} = ['-- 仿射等价类：规模超限，已跳过' ...
            '（增大 opts.MaxOps 可强制计算）。--'];
    else
        L{end+1} = sprintf(['-- 仿射等价类（B ~ u*B+t，gcd(u,%d)=1；' ...
            '共 %d 类，代表如下）--'], Mp, cls.affine.count);
        L = emit_reps(L, cls.affine.reps, opts.ShowReps);
    end

    % ---- 砖块分拆等价类
    L{end+1} = '';
    if isfield(cls, 'partitionSkipped') && cls.partitionSkipped
        L{end+1} = ['-- 砖块分拆等价类：铺法总数过大，已跳过' ...
            '（增大 opts.MaxOps 可强制计算）。--'];
    else
        L{end+1} = sprintf(['-- 砖块分拆等价类（equivalent as tile partitions，' ...
            '{A+b : b∈B} 相同；对全部 %d 个铺法计算；共 %d 类，代表如下）--'], ...
            out.numTilings, cls.partition.count);
        L = emit_reps(L, cls.partition.reps, opts.ShowReps);
    end

    L{end+1} = '';

    S.periods{t} = struct('M', Mp, 'out', out, 'classify', cls);
end

txt = strjoin(L, newline);
end

% =========================================================================
function L = emit_reps(L, reps, showMax, perReps)
% 列出类代表；perReps 可选（与 reps 一一对应的基本周期标注）
if nargin < 4, perReps = {}; end
nr = numel(reps);
show = min(nr, showMax);
for r = 1:show
    tag = '';
    if ~isempty(perReps) && ~isempty(perReps{r})
        tag = sprintf('   （基本周期 %d）', perReps{r});
    end
    L{end+1} = sprintf('  %3d)  B = %s%s', r, fmt_set(reps{r}), tag); %#ok<AGROW>
end
if nr > show
    L{end+1} = sprintf('  ... 其余 %d 个代表略（结果结构体中完整保存）。', nr - show);
end
if nr == 0
    L{end+1} = '  （无）';
end
end

% =========================================================================
function v = layerPeriods(out)
% 每个代表所属层的基本周期（ = P * table(t,1) ）；repsByLayer 按层序拼接
% 恰为 out.representatives 的顺序
v = cell(1, numel(out.representatives));
pos = 0;
for t = 1:size(out.table, 1)
    Pd = out.table(t, 2);
    rl = out.repsByLayer{t};
    for r = 1:numel(rl)
        pos = pos + 1;
        if pos <= numel(v)
            v{pos} = Pd;
        end
    end
end
end

% =========================================================================
function s = class_str(c)
if isnan(c), s = '（跳过）'; else, s = sprintf('%d', c); end
end

% =========================================================================
function d = all_divisors(n)
n = double(n);
d = 1;
f = reshape(double(factor(n)), 1, []);
for i = 1:numel(f)
    d = unique([d, d * f(i)]);
end
d = sort(reshape(d, 1, []));
end

% =========================================================================
function s = numjoin(v)
v = v(:)';
s = sprintf('%d, ', v);
if ~isempty(s), s = s(1:end-2); end
end

% =========================================================================
function s = fmt_set(B)
if isempty(B), s = '[]'; return; end
s = ['[' strtrim(sprintf('%d ', B)) ']'];
end
