function [A2, M2, info] = reduce_tiling_period(A, M, opts)
%REDUCE_TILING_PERIOD  把 (A, M) 约化到"基周期"：M → M/|H|
%
%   [A2, M2, info] = reduce_tiling_period(A, M)
%   [A2, M2, info] = reduce_tiling_period(A, M, opts)
%
% 一、A 是子群陪集之并（含 A 本身为子群）
%   设 Stab(A) = { h ∈ Z_M : A + h = A }，它是 Z_M 的子群，阶为 m。
%   A 是某子群若干陪集之并  <=>  m > 1。取 d = M/m（此时 Stab(A) = {0,d,2d,…}，
%   即 d 的所有倍数），则 A = π^{-1}(Ā)（π: Z_M → Z_d 为自然投影），于是
%
%        A ⊕ B = Z_M   <=>   Ā ⊕ B̄ = Z_d ,  Ā = A mod d , B̄ = B mod d
%
%   且每个 Z_d 的铺法 B̄ 恰好对应 m^{|B|} 个 Z_M 的铺法（B 的每个元素可自由
%   加上 H 中元素），所以  F_A(Z_M) = m^{|B|} · F_A(Z_d)。
%   直观上：彩虹图沿圆周每 d 个扇区完全重复 m 次。
%
% 二、A 无周期，但 M 仍不是基周期（存在陪集并的补集）
%   若存在真因子 d | M 使得 A mod d 仍是 |A| 个不同元素（A 在 Z_d 中是横截），
%   且 Ā = A mod d 能密铺 Z_d，则把 (A, M) 换成 (Ā, d)。此时
%        { Z_M 中满足 B 为陪集并的补集 }  ←→  { Z_d 的补集 }
%   是一一对应（B ↔ B̄ = B mod d）；Z_M 中若还有其它补集，它们不属于这一族。
%
% 算法：先做一（结构判定，无需枚举），再做二（对 M 的真因子升序扫描，取最小 d），
%       反复直到两者都不再成立；终止条件另加：|A| = 1 或 |A| = M（图形已平凡）。
%
% 选项
%   opts.MaxSteps  最多约化步数（默认 30）
%   opts.MaxEnum   内部判定"能否密铺 Z_d"时的枚举上限（默认 2e5）
%
% 输出
%   A2, M2 : 约化后的密铺集与模数（A2 ⊂ [0, M2)，|A2| 仍为 p 的幂）
%   info   : .reduced .A0 .M0 .A .M .factor .steps{kinds,orders}
%            .nTilingsFactor（仅"A 型"约化累计的 F 倍数）
%            .checked（是否做过检查，由调用者置位）

if nargin < 3 || isempty(opts), opts = struct(); end
dflt = struct('MaxSteps', 30, 'MaxEnum', 2e5);
fn = fieldnames(dflt);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = dflt.(fn{i}); end
end

A2 = sort(unique(mod(double(A(:)'), double(M))));
M2 = double(M);

info = struct('reduced', false, 'A0', A2, 'M0', M2, 'A', A2, 'M', M2, ...
              'factor', 1, 'steps', {{}}, 'kinds', {{}}, 'orders', [], ...
              'nTilingsFactor', 1, 'dmods', M2);

for it = 1:opts.MaxSteps
    if numel(A2) <= 1 || numel(A2) >= M2 || M2 <= 1
        break                                   % 图形已平凡，不再约化
    end

    % ---------- 一、A 自身是陪集并 ----------
    [mA, ~] = period_order(A2, M2);
    if mA > 1
        d = M2 / mA;
        if numel(A2) == mA
            extra = '（A 本身就是这个子群）';
        else
            extra = sprintf('（A 是 H 的 %d 个陪集之并）', numel(A2) / mA);
        end
        info.steps{end+1} = sprintf( ...
            ['A 以阶 %d 的子群 H = {%s} 为周期%s ' ...
             '→ 商掉后 M: %d → %d（彩虹图沿圆周重复 %d 次）'], ...
            mA, str_set(0:d:M2-1, 6), extra, M2, d, mA);
        info.kinds{end+1}  = 'A';
        info.orders(end+1) = mA;
        info.nTilingsFactor = info.nTilingsFactor * mA ^ (M2 / numel(A2));
        A2 = sort(unique(mod(A2, d)));
        M2 = d;
        info.dmods(end+1) = M2;
        continue
    end

    % ---------- 二、M 不是基周期：存在更小模数上可密铺的横截 ----------
    ds = my_divisors(M2);
    hit = false;
    for i = 1:numel(ds)
        d = ds(i);
        if d >= M2, break; end
        Ad = sort(unique(mod(A2, d)));
        if numel(Ad) ~= numel(A2)
            continue                            % A mod d 不是横截，跳过
        end
        if tiles_in(Ad, d, opts.MaxEnum)
            info.steps{end+1} = sprintf( ...
                ['A 无周期，但 A mod %d 仍是 %d 个不同元素，且能密铺 Z_%d ' ...
                 '（对应的补集是阶 %d 子群的陪集之并）→ M: %d → %d'], ...
                d, numel(A2), d, M2 / d, M2, d);
            info.kinds{end+1}  = 'B';
            info.orders(end+1) = M2 / d;
            A2 = Ad;
            M2 = d;
            info.dmods(end+1) = M2;
            hit = true;
            break
        end
    end
    if ~hit, break; end
end

info.A = A2;
info.M = M2;
info.factor = info.M0 / M2;
info.reduced = (M2 < info.M0);
end

% -------------------------------------------------------------------------
function tf = tiles_in(X, d, maxEnum)
% X 能否密铺 Z_d（调用主引擎；不绘图、不打印）
tf = false;
if numel(X) > d, return; end
try
    o = tile_all_tilings(X, d, struct('Report', false, 'MaxEnum', maxEnum));
    tf = o.isTile && o.numClasses >= 1;
catch
    tf = false;
end
end

% -------------------------------------------------------------------------
function [m, h] = period_order(X, M)
% X 在 Z_M 中稳定子 Stab(X) = {h : X+h = X} 的阶 m 与生成元 h = M/m。
% 稳定子是子群，而 Z_M 的子群都由 M/m（m | M）生成，故只需对 m 的各因子检验。
X = sort(unique(mod(double(X(:)'), M)));
m = 1;  h = 0;
if isempty(X), return; end
if numel(X) >= M, m = M; h = 1; return; end
ds = my_divisors(M);
for i = 1:numel(ds)
    mi = ds(i);
    if mi <= m, continue; end
    if isequal(sort(mod(X + M / mi, M)), X)
        m = mi;
        h = M / mi;
    end
end
end

% -------------------------------------------------------------------------
function d = my_divisors(n)
d = 1;
f = reshape(double(factor(n)), 1, []);
for i = 1:numel(f)
    d = unique([d, d * f(i)]);
end
d = sort(d);
end

% -------------------------------------------------------------------------
function s = str_set(v, maxItems)
v = double(v(:)');
if nargin < 2, maxItems = 10; end
if numel(v) > maxItems
    s = sprintf('%s, ...（共 %d 个）', strtrim(sprintf('%d, ', v(1:maxItems))), numel(v));
else
    s = strtrim(sprintf('%d, ', v));
end
if ~isempty(s) && s(end) == ',', s(end) = []; end
end
