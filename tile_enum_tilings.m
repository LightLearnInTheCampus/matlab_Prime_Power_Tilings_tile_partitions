function out = tile_enum_tilings(A, M, opts)
%TILE_ENUM_TILINGS  枚举素数幂基数密铺集 A 在 Z_M 中的全部铺法（按循环平移分类）
%   26-11-17 版本：不带任何绘图；报告文本不含"实际加载"路径。
%
%   out = tile_enum_tilings(A)           % M 缺省取万能（最大基本）周期 P*D
%   out = tile_enum_tilings(A, M)        % 指定模数 M（必须满足 P | M）
%   out = tile_enum_tilings(A, M, opts)  % 附加选项
%
% 输入
%   A    : 整数向量（允许负元素，内部平移规范化），|A| = p^alpha，p 为素数
%   M    : 模数；为空 [] 时自动取 P*D
%   opts : 结构体（可缺省）
%          .ReducePeriod 是否先判定 A（或补集 B）是否为子群陪集之并，
%                        并据此把 M 约化到基周期后再求解        (默认 false)
%          .MaxEnum      单层枚举上限（超出则只计数）(默认 2e5)
%          .Report       是否打印报告              (默认 true)
%          .ShowReps     报告中最多个列出的类代表   (默认 40)
%
% 输出 out
%   .isTile          A 是否为整数密铺集
%   .p .k .nu        p，k = log_p|A|，素幂层 nu_1<...<nu_k
%   .Delta           fiber gcd（Delta_1 | ... | Delta_k = D）
%   .P .D .PD        P = p^{nu_k}，D = Delta_k，万能周期 PD
%   .canonicalTile   典范（补集等价）密铺集 A^*
%   .M .dmax         M 与 gcd(M/P, D)
%   .table           [d, Pd, F_A(Pd), E_A(Pd), 平移类数]
%   .numTilings      Z_M 中的铺法总数
%   .numClasses      循环平移等价类数
%   .representatives 每个平移类的一个代表（模 M，相对原始输入 A）
%   .repsByLayer{t}  第 t 层（基本周期 P*table(t,1)）的平移类代表（已膨胀到 Z_M）
%   .complete        是否完成了完整枚举（受 MaxEnum 限制）
%   .report          报告文本
%
% 理论依据（论文 "Fundamental periods and complements of prime-power
% integer tiles"）：
%   Theorem 3.1   基本周期恰为 {Pd : d | D}，且 P|M 时 Z_M 才有铺法
%   Theorem 4.2   截面塔定理：所有补集由逐层截面函数 alpha_i: C_i -> F_p
%                 唯一生成，F_A(Pd) = p^{sum_i p^{nu_i-i} gcd(d,Delta_i)}
%   Corollary 4.4 E_A(Pd) 由 Moebius 反演给出，平移类数 = E_A(Pd)/(Pd)
%   Prop. 4.9     精确基本周期判定：B + Pd/q ~= B 对每个素因子 q | d
%
% 例：
%   out = tile_enum_tilings([0 3 10 13]);                % P=4, D=5, M=20
%   out = tile_enum_tilings([0 3 6], 9);
%
% 见 tile_enum_tilings_test.m 中的论文算例校验。

% ---------------------------------------------------------------- 0. 选项
TILE_VER = '2026-11-17a';   % 版本戳：26-11-17 版本（无绘图）
if nargin < 2, M = []; end
if nargin < 3 || isempty(opts), opts = struct(); end
dflt = struct('MaxEnum', 2e5, 'Report', true, 'ShowReps', 40, ...
              'ReducePeriod', false);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = dflt.(fn{i}); end
end

% ------------------------------------------------- 1. 输入检查与平移规范化
if ~isnumeric(A) || isempty(A)
    error('tile_enum_tilings:input', 'A 必须是非空整数向量。');
end
A = double(A(:)');
if any(~isfinite(A)) || any(A ~= round(A))
    error('tile_enum_tilings:input', 'A 的元素必须为整数。');
end
A = sort(A);
if numel(unique(A)) ~= numel(A)
    error('tile_enum_tilings:input', 'A 的元素必须两两不同。');
end
n    = numel(A);
amin = min(A);
A0   = A - amin;          % 规范形：min A0 = 0（平移不改变铺法族）

out = struct('input', A, 'A0', A0, 'n', n, 'isTile', false, ...
             'version', TILE_VER, ...
             'p', NaN, 'k', NaN, 'nu', [], 'Delta', [], 'P', NaN, ...
             'D', NaN, 'PD', NaN, 'canonicalTile', [], ...
             'M', [], 'dmax', [], 'table', zeros(0,5), ...
             'numTilings', 0, 'numClasses', 0, 'representatives', {{}}, ...
             'repsByLayer', {{}}, ...
             'complete', true, 'report', '', ...
             'reduce', struct('reduced', false, 'checked', false));

% ------------- 1b. 可选：先判定 A（或补集 B）是否为陪集并，约化到基周期
if opts.ReducePeriod && ~isempty(M)
    [A2, M2, rinfo] = reduce_tiling_period(A0, M, struct('MaxEnum', opts.MaxEnum));
    rinfo.checked = true;
    out.reduce = rinfo;
    if rinfo.reduced
        out.M0 = M;                       % 约化前的模数
        A0 = A2;  A = A2;  amin = 0;  n = numel(A0);  M = M2;
        out.A0 = A0;  out.n = n;
    end
end
diam = max(A0);

% ------------------------------------------------------ 2. n = 1 平凡情形
if n == 1
    if isempty(M), M = 1; end
    validate_M(M);
    out.isTile = true;  out.p = 1;  out.k = 1;  out.nu = 0;
    out.Delta = 1;  out.P = 1;  out.D = 1;  out.PD = 1;
    out.canonicalTile = A0;
    out.M = M;  out.dmax = 1;
    out.table = [1, 1, 1, 1, 1];
    out.numTilings = 1;  out.numClasses = 1;
    out.representatives = {0:M-1};
    out.repsByLayer = {{0:M-1}};
    out.report = report_text(out, opts);
    if opts.Report, fprintf('%s\n', out.report); end
    return
end

% --------------------------------------------------- 3. 单素数幂分解 |A|=p^k
fac = factor(n);
p   = fac(1);
if any(fac ~= p)
    error('tile_enum_tilings:notPrimePower', ...
        '|A| = %d 不是单素数幂 p^alpha；本程序只处理 |A| = p^alpha。', n);
end
k = numel(fac);

% ---------------------------------------- 4. 素幂层：Phi_{p^j} | A(x) 的 j
mask = zeros(1, diam + 1);
mask(A0 + 1) = 1;                       % 规范化 mask 多项式（系数 0/1）
lev = zeros(1, 0);
q   = 1;
while (p - 1) * q <= diam
    phi = zeros(1, (p - 1) * q + 1);    % Phi_{p^j} = 1+x^q+...+x^{(p-1)q}
    phi(1:q:end) = 1;                   % 系数 1 落在指数 0,q,2q,...,(p-1)q
    if poly_divisible_exact(mask, phi)
        lev(end + 1) = p * q;           %#ok<AGROW>
    end
    q = q * p;
end

if numel(lev) ~= k
    % 素数幂判据：Phi_{p^j} 因子个数必须恰为 k = log_p|A|
    out.report = sprintf(['A = %s 不是整数密铺集：\n' ...
        '  |A| = %d = %d^%d，但 mask 中被 Phi_{p^j} 整除的层只有 %d 个（应为 %d 个）。\n' ...
        '  按 (T1) 判据，A 不能密铺 Z（因而也不能密铺任何 Z_M）。'], ...
        fmt_set(A), n, p, k, numel(lev), k);
    if opts.Report, fprintf('%s\n', out.report); end
    return
end

nu = round(log(lev) / log(p));          % nu_i = log_p(lev_i)
nuEnd = nu(end);
P     = p ^ nuEnd;

% -------------------------------------------------------- 5. fiber gcd Delta_i
Delta = zeros(1, k);
for i = 1:k
    Qi  = p ^ (nu(i) - 1);              % Q_i = p^{nu_i-1}
    cls = mod(A0, Qi);
    uc  = unique(cls);
    g   = 0;
    for c = uc(:)'
        idx = find(cls == c);
        if numel(idx) < 2, continue; end
        dd = A0(idx) - A0(idx(1));      % 该类内相对第一个元素的差
        g  = gcd(g, gcd_vec(dd) / Qi);  % 全部差都是 Q_i 的倍数
        if g == 1, break; end
    end
    Delta(i) = g;
end
D = Delta(k);

% -------------------------------------------------- 6. 典范（补集等价）密铺集
gstep = p .^ (nu - 1) .* Delta;         % g_i = p^{nu_i-1} Delta_i
Astar = 0;
for i = 1:k
    Ai    = (0:p-1) * gstep(i);
    Astar = unique(reshape(Astar(:) + Ai(:).', 1, []));
end
Astar = sort(Astar);

out.isTile = true;  out.p = p;  out.k = k;  out.nu = nu;
out.Delta = Delta;  out.P = P;  out.D = D;  out.PD = P * D;
out.canonicalTile = Astar;

% ------------------------------------------------------------- 7. 模数 M
if isempty(M)
    M = P * D;                          % 万能模数（最大基本周期）
end
validate_M(M);
if mod(M, P) ~= 0
    out.M = M;  out.report = sprintf([ ...
        '模数 M = %d 不是 P = p^{nu_k} = %d 的倍数，\n' ...
        '由论文 Lemma 3.2 / Theorem 3.1，P 整除每个补集的每个周期，故 Z_%d 中无铺法。'], ...
        M, P, M);
    if opts.Report, fprintf('%s\n', out.report); end
    return
end
s    = M / P;
dmax = gcd(s, D);                       % gcd(M/P, D)：所有可达的基本周期

out.M = M;  out.dmax = dmax;

% --------------------------------------------------- 8. 逐层（按基本周期）枚举
ds   = divisors_of(dmax);               % 所有 d | dmax，升序
nS   = numel(ds);
tab  = zeros(nS, 5);
repsM = {};                             % 模 M 下的类代表（已平移回 A 的坐标）
repsByLayer = cell(1, nS);              % 第 t 层（基本周期 P*ds(t)）的类代表
cntEnum = 0;
for t = 1:nS
    d = ds(t);
    N = P * d;                          % 该层模数
    Fd = p ^ G_A(d, p, nu, Delta);      % F_A(Pd)：以 Pd 为周期的补集个数
    Ed = 0;                             % E_A(Pd)：基本周期恰为 Pd 的个数
    layerReps = {};                     % 本层（膨胀到 Z_M 后）的类代表
    for e = divisors_of(d)
        Ed = Ed + mobius_of(d / e) * p ^ G_A(e, p, nu, Delta);
    end
    tab(t, 1:3) = [d, N, Fd];
    tab(t, 4)   = Ed;

    if Fd == 0 || Ed == 0
        tab(t, 5) = 0;
        continue
    end

    if cntEnum + Fd > opts.MaxEnum
        % 输出规模过大：只给计数（论文 §6.2 的建议：计数先行）
        tab(t, 5) = Ed / N;
        out.complete = false;
        warning('tile_enum_tilings:tooMany', ...
            ['模数 %d 的补集数 %d 超过 MaxEnum = %d，该层只给出计数、未枚举。\n' ...
             '如需完整列表，请增大 opts.MaxEnum（内存与时间随规模增长）。'], ...
            N, Fd, opts.MaxEnum);
        continue
    end

    try
        Bmat = tower_sets(d, p, nu, Delta, k, opts.MaxEnum);  % 该层全部补集
        if size(Bmat, 1) ~= Fd
            error('tile_enum_tilings:internal', ...
                '内部校验失败：枚举数 %d != 公式 F_A(%d) = %d。', size(Bmat,1), N, Fd);
        end
        Bex  = Bmat(exact_period_mask(Bmat, N), :);          % 基本周期恰为 N
        if size(Bex, 1) ~= Ed
            error('tile_enum_tilings:internal', ...
                '内部校验失败：基本周期层大小 %d != 公式 E_A(%d) = %d。', size(Bex,1), N, Ed);
        end
        [repsN, ~] = classes_by_translation(Bex, N);         % 平移等价类
        tab(t, 5)  = size(repsN, 1);

        % 膨胀到 Z_M，并平移回原始 A 的坐标：A0 + B = Z_M  <=>  A + (B-amin) = Z_M
        mult = M / N;
        srep = size(repsN, 2);
        for r = 1:size(repsN, 1)
            Bm = zeros(1, srep * mult);           % 显式拼装，不用隐式扩展
            for o = 0:mult-1
                Bm(o * srep + (1:srep)) = repsN(r, :) + N * o;
            end
            Bm = sort(unique(mod(Bm, M)));
            Bm = mod(Bm - amin, M);
            % 规范化：平移使 0 ∈ B（B-min(B) 仍是同一平移类的合法铺法）
            Bm = mod(Bm - min(Bm), M);
            repsM{end + 1} = sort(Bm);                         %#ok<AGROW>
            layerReps{end + 1} = sort(Bm);                     %#ok<AGROW>
        end
        repsByLayer{t} = layerReps;
        cntEnum = cntEnum + Fd;
    catch ME
        fprintf(2, '\n[调试信息] 第 d = %d 层（N = %d）出错：%s\n', d, N, ME.message);
        fprintf(2, '  p = %g, k = %d, nu = %s, Delta = %s, P = %g, D = %g, dmax = %g\n', ...
            p, k, mat2str(nu), mat2str(Delta), P, D, dmax);
        fprintf(2, '  F_A(Pd) = %g, E_A(Pd) = %g, M = %g\n', Fd, Ed, M);
        if exist('Bmat', 'var') %#ok<EXIST>
            fprintf(2, '  size(Bmat) = %s, class(Bmat) = %s\n', ...
                mat2str(size(Bmat)), class(Bmat));
        end
        rethrow(ME)
    end
end

% ------------------------------------------------ 8.5 层内字典序排序
% 每个基本周期层的代表按字典序排列（sortrows 逐行比较），并按层序
% 重建 representatives 与 repsByLayer，保证输出与层归属一致。
for t = 1:nS
    rl = repsByLayer{t};
    if isempty(rl), continue; end
    nL  = numel(rl);
    mat = zeros(nL, numel(rl{1}));
    for r = 1:nL
        mat(r, :) = rl{r};
    end
    mat = sortrows(mat);                 % 字典序（逐列比较，元素均在 Z_M 内）
    for r = 1:nL
        rl{r} = mat(r, :);
    end
    repsByLayer{t} = rl;
end
repsM = cell(1, 0);
for t = 1:nS
    for r = 1:numel(repsByLayer{t})
        repsM{end + 1} = repsByLayer{t}{r};              %#ok<AGROW>
    end
end

out.table          = tab;
out.numTilings     = p ^ G_A(M / P, p, nu, Delta);   % F_A(M) 的公式值
out.numClasses     = sum(tab(:, 5));
out.representatives = repsM;
out.repsByLayer    = repsByLayer;

% 交叉校验：F_A(M) 应等于各层 E_A(Pd) 之和
if sum(tab(:,4)) ~= out.numTilings
    error('tile_enum_tilings:internal', ...
        '内部校验失败：sum E_A(Pd) = %d != F_A(M) = %d。', sum(tab(:,4)), out.numTilings);
end

% ------------------------------------------------------------- 9. 报告
out.report = report_text(out, opts);
if opts.Report, fprintf('%s\n', out.report); end
end
% ================================================================ 主函数结束


% -------------------------------------------------------------------------
function validate_M(M)
if ~isscalar(M) || ~isfinite(M) || M ~= round(M) || M < 1
    error('tile_enum_tilings:badM', '模数 M 必须是正整数。');
end
end

% -------------------------------------------------------------------------
function G = G_A(d, p, nu, Delta)
% 论文 (4.2)：G_A(d) = sum_i p^{nu_i - i} gcd(d, Delta_i)
G = 0;
for i = 1:numel(nu)
    G = G + p ^ (nu(i) - i) * gcd(d, Delta(i));
end
end

% -------------------------------------------------------------------------
function tf = poly_divisible_exact(a, b)
% 整数多项式整除判定：b | a ?
% 约定 a、b 都是升幂系数向量（下标 k 对应 x^{k-1}），且 b 的常数项 b(1) = 1。
% 因 b(1) = 1，可从低次向高次逐项消去，余项全零即整除；全程整数运算，无容差。
la = numel(a); lb = numel(b);
if lb > la
    tf = false; return
end
r = a;
for k = 1:(la - lb + 1)
    c = r(k);
    if c ~= 0
        r(k:k + lb - 1) = r(k:k + lb - 1) - c * b;   % 消去 x^{k-1} 项
    end
end
tf = all(r == 0);
end

% -------------------------------------------------------------------------
function g = gcd_vec(v)
g = 0;
for i = 1:numel(v)
    g = gcd(g, abs(v(i)));
    if g == 1, return; end
end
end

% -------------------------------------------------------------------------
function d = divisors_of(n)
% n 的全部正因数（升序行向量）
% 注意：factor(n) 返回"行"向量，若写成 for f = factor(n)'，for 是按列遍历，
% 对 2x1 列只有一列，f 会拿到整个列向量 —— 必须用下标显式取标量。
n = double(n);
n = n(1);                                   % 强制标量
d = 1;
f = reshape(double(factor(n)), 1, []);      % 素因子统一成行向量
for i = 1:numel(f)
    d = unique([d, d * f(i)]);              %#ok<AGROW>
end
d = sort(reshape(d, 1, []));                % 返回行向量
end

% -------------------------------------------------------------------------
function m = mobius_of(n)
if n == 1, m = 1; return; end
f = factor(n);
if numel(unique(f)) < numel(f), m = 0; return; end
m = (-1) ^ numel(f);
end

% -------------------------------------------------------------------------
function Bmat = tower_sets(d, p, nu, Delta, k, maxEnum)
% 论文 Theorem 4.2 的截面塔：生成模数 N = p^{nu_k} d 下、以 N 为周期的
% 全部补集，返回矩阵（每行一个升序集合）。第 i 层：
%   d_i = gcd(d, Delta_i),  N_i = p^{nu_i} d_i,  L_i = N_i / p
%   B_i = { c + L_i * alpha_i(c) : c in C_i },  C_i = B_{i-1} 在 Z_{L_i} 的原像
Bmat  = 0;          % B_0 = {0} \subset Z_1
Nprev = 1;
for i = 1:k
    di = gcd(d, Delta(i));
    Ni = p ^ nu(i) * di;
    Li = Ni / p;
    if mod(Li, Nprev) ~= 0
        error('tile_enum_tilings:internal', '内部错误：N_{i-1} 不整除 L_i。');
    end
    szB   = size(Bmat, 2);
    nC    = szB * (Li / Nprev);          % |C_i| = p^{nu_i-i} d_i
    total = p ^ nC;                      % 该层截面函数个数 p^{|C_i|}
    if size(Bmat, 1) * total > maxEnum
        error('tile_enum_tilings:tooMany', ...
            '枚举规模 %d 超过 MaxEnum = %d。', size(Bmat,1) * total, maxEnum);
    end
    at = digits_table(total, p, nC);     % 所有 alpha 的 p 进制数字表
    nb = cell(size(Bmat, 1), 1);
    for r = 1:size(Bmat, 1)
        C = inflate_set(Bmat(r, :), Nprev, Li);
        if numel(C) ~= nC
            error('tile_enum_tilings:internal', '内部错误：|C_i| 不符。');
        end
        tmp  = mod(repmat(C, total, 1) + Li * at, Ni);
        nb{r} = sort(tmp, 2);            % 规范代表（升序）
    end
    Bmat  = vertcat(nb{:});
    Nprev = Ni;
end
end

% -------------------------------------------------------------------------
function C = inflate_set(B, Nprev, L)
% B \subset Z_{Nprev} 在 Z_L 中的原像（Nprev | L）
m  = L / Nprev;
B  = B(:)';
nb = numel(B);
C  = zeros(1, nb * m);                  % 显式拼装，避免依赖隐式扩展
for j = 0:m-1
    C(j * nb + (1:nb)) = B + Nprev * j;
end
C = unique(mod(C, L));
end

% -------------------------------------------------------------------------
function D = digits_table(total, p, width)
% 返回 total x width 矩阵，第 r 行是 (r-1) 的 p 进制数字（低位在前）
if width == 0, D = zeros(total, 0); return; end
v = (0:total-1)';
D = zeros(total, width);
for j = 1:width
    D(:, j) = mod(floor(v / p ^ (j-1)), p);
end
end

% -------------------------------------------------------------------------
function msk = exact_period_mask(Bmat, N)
% Prop. 4.9：B 的基本周期为 N  <=>  对每个素因子 q | N，B + N/q ~= B
% 逐行、逐素因子计算，不用任何隐式扩展，任何形状/类型都安全。
F   = size(Bmat, 1);
msk = true(F, 1);
N   = double(N(1));
if N <= 1, return; end
qlist = unique(factor(N));
qlist = reshape(double(qlist), 1, []);      % 强制成行向量，逐个取标量
for j = 1:numel(qlist)
    K = N / qlist(j);                       % 标量平移量
    for r = 1:F
        if ~msk(r), continue; end
        row  = Bmat(r, :);                  % 1 x s
        srow = sort(mod(row + K, N));       % 1 x s，平移后仍升序
        if all(srow == row)                 % B + K == B，周期不精确
            msk(r) = false;
        end
    end
end
end

% -------------------------------------------------------------------------
function [reps, orbit] = classes_by_translation(Bmat, N)
% 按循环平移 B ~ B + t 分类，返回代表（矩阵行）与各轨道大小
F = size(Bmat, 1);
if F == 0
    reps = zeros(0, size(Bmat, 2)); orbit = []; return
end
if N <= 63
    pow2v = zeros(N, 1, 'uint64');
    for b = 0:N-1
        pow2v(b + 1) = bitshift(uint64(1), b);   % 2^b，标量移位，避免类型混用
    end
    code = zeros(F, 1, 'uint64');
    for j = 1:size(Bmat, 2)
        code = bitor(code, pow2v(Bmat(:, j) + 1));
    end
    cmin = min_rotation_code(code, N);
    [~, ia, ic] = unique(cmin);
    reps  = Bmat(ia, :);
    orbit = accumarray(ic, 1);
    if any(orbit ~= N)
        error('tile_enum_tilings:internal', ...
            '内部错误：平移轨道大小不是 N = %d（周期不精确？）。', N);
    end
else
    % N > 63：退化为逐集合比较（极慢，一般不会出现）
    reps = zeros(0, size(Bmat, 2)); orbit = [];
    keys = cell(F, 1);
    for r = 1:F
        keys{r} = rotation_key(Bmat(r, :), N);
    end
    [~, ia, ic] = unique(keys);
    reps  = Bmat(ia, :);
    orbit = accumarray(ic, 1);
    if any(orbit ~= N)
        error('tile_enum_tilings:internal', '内部错误：平移轨道大小不是 N。');
    end
end
end

% -------------------------------------------------------------------------
function c = min_rotation_code(c, M)
% 长度为 M 的循环二进制串的最小旋转（code 用 uint64 位表示，M <= 63）
if M <= 1, return; end
mask = bitshift(uint64(1), M) - uint64(1);
top  = bitshift(uint64(1), M - 1);
cur  = c;
for t = 1:M-1
    carry = bitand(cur, top) ~= 0;
    cur   = bitand(bitshift(cur, 1), mask);
    if any(carry)
        cur(carry) = bitor(cur(carry), uint64(1));
    end
    c = min(c, cur);
end
end

% -------------------------------------------------------------------------
function k = rotation_key(B, N)
% 平移不变量：轨道上所有集合"排序后元组"的字典序最小者（与 canon_key 同一约定）
B = sort(mod(B(:)', N));
t = numel(B);
if t == 0, k = ''; return; end
R = zeros(N, t);
for s = 1:N
    R(s, :) = sort(mod(B + s - 1, N));
end
R = sortrows(R, 1:t);
k = mat2str(R(1, :));
end

% -------------------------------------------------------------------------
function s = fmt_set(B)
if isempty(B), s = '[]'; return; end
s = ['[' strtrim(sprintf('%d ', B)) ']'];
end

% -------------------------------------------------------------------------
function txt = report_text(out, opts)
L = {};
L{end+1} = '================ 密铺集 A 的算术数据 ================';
% ---- 基周期约化说明（若启用）----
if isfield(out, 'reduce') && isfield(out.reduce, 'steps') && ~isempty(out.reduce.steps)
    L{end+1} = '-------- 基周期约化：A 或补集是子群陪集之并 --------';
    for i = 1:numel(out.reduce.steps)
        L{end+1} = sprintf('  第 %d 步：%s', i, out.reduce.steps{i}); %#ok<AGROW>
    end
    L{end+1} = sprintf('  M 缩小 %d 倍：Z_%d → Z_%d ；约化后 A* = %s（|A*| = %d）', ...
        out.reduce.factor, out.reduce.M0, out.reduce.M, ...
        fmt_set(out.reduce.A), numel(out.reduce.A));
    if isfield(out.reduce, 'kinds') && any(strcmp(out.reduce.kinds, 'B'))
        L{end+1} = sprintf(['  注：补集为陪集并，Z_%d 中"补集是子群陪集并"的那一族铺法' ...
            '与 Z_%d 的铺法一一对应（可由后者完全提升得到）；' ...
            'Z_%d 若还有其它补集，不在本次约化范围内。以下数据按 Z_%d 给出。'], ...
            out.reduce.M0, out.reduce.M, out.reduce.M0, out.reduce.M);
    else
        L{end+1} = sprintf(['  注：A 为陪集并，Z_%d 的铺法数 = %d × (Z_%d 的铺法数)；' ...
            '以下数据均按约化后的 Z_%d 给出。'], ...
            out.reduce.M0, out.reduce.nTilingsFactor, out.reduce.M, out.reduce.M);
    end
    L{end+1} = '-----------------------------------------------------';
elseif isfield(out, 'reduce') && isfield(out.reduce, 'checked') && out.reduce.checked
    L{end+1} = '已检查基周期：A 与任一补集都不是子群陪集之并，M 已是基周期。';
end
if out.n == 1
    if isfield(out, 'reduce') && out.reduce.reduced
        L{end+1} = sprintf('约化后 A* = %s（|A*| = 1），Z_%d 中唯一铺法 B = [0..%d]，1 个平移类。', ...
            fmt_set(out.A0), out.M, out.M - 1);
    else
        L{end+1} = sprintf('A = %s   （|A| = 1，平凡情形）', fmt_set(out.input));
        L{end+1} = sprintf('模数 M = %d：唯一铺法 B = [0..%d]，1 个平移类。', out.M, out.M-1);
    end
    txt = strjoin(L, newline); return
end
if ~out.isTile
    L{end+1} = sprintf('A = %s', fmt_set(out.input));
    L{end+1} = '结论：A 不是整数密铺集，无任何铺法。';
    txt = strjoin(L, newline); return
end
if isfield(out, 'reduce') && isfield(out.reduce, 'reduced') && out.reduce.reduced
    L{end+1} = sprintf('A = %s   （原始输入；约化后 A* = %s）', ...
        fmt_set(out.input), fmt_set(out.A0));
else
    L{end+1} = sprintf('A = %s   （规范化 A0 = %s，min = 0）', ...
        fmt_set(out.input), fmt_set(out.A0));
end
L{end+1} = sprintf('|A| = %d = %d^%d        p = %d,  alpha = k = %d', ...
    out.n, out.p, out.k, out.p, out.k);
L{end+1} = sprintf('素幂层  nu_i = %s        P = p^{nu_k} = %d', ...
    fmt_set(out.nu), out.P);
L{end+1} = sprintf('fiber gcd  Delta_i = %s   D = Delta_k = %d', ...
    fmt_set(out.Delta), out.D);
L{end+1} = sprintf('基本周期集合  {P*d : d | D} 中最大者  PD = %d', out.PD);
L{end+1} = sprintf('典范等价密铺集  A* = %s', fmt_set(out.canonicalTile));
L{end+1} = '-----------------------------------------------------';
L{end+1} = sprintf('模数 M = %d   （M/P = %d,  gcd(M/P, D) = %d）', ...
    out.M, out.M / out.P, out.dmax);
L{end+1} = '';
L{end+1} = '   d      Pd     F_A(Pd)     E_A(Pd)   平移类数';
L{end+1} = '  ------------------------------------------------';
for t = 1:size(out.table, 1)
    L{end+1} = sprintf('  %2d  %6d  %10d  %10d  %8d', out.table(t, :)); %#ok<AGROW>
end
L{end+1} = '  ------------------------------------------------';
L{end+1} = sprintf('合计：Z_%d 中铺法总数 F_A(M) = %d', out.M, out.numTilings);
L{end+1} = sprintf('      循环平移等价类数        = %d', out.numClasses);
if ~out.complete
    L{end+1} = '注意：部分层的枚举被 MaxEnum 限制跳过（计数仍然精确）。';
end
L{end+1} = '-----------------------------------------------------';
L{end+1} = sprintf('每个平移类等价于一个环状彩虹图；下面列出各类代表（模 %d）：', out.M);
nr = numel(out.representatives);
show = min(nr, opts.ShowReps);
for r = 1:show
    L{end+1} = sprintf('  %3d)  B = %s', r, fmt_set(out.representatives{r})); %#ok<AGROW>
end
if nr > show
    L{end+1} = sprintf('  ... 其余 %d 个类代表见 out.representatives', nr - show);
end
txt = strjoin(L, newline);
end
