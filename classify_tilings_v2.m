function out = classify_tilings_v2(A, M, reps, opts)
%CLASSIFY_TILINGS_V2  对 Z_M 中密铺集 A 的各平移类代表再做两种 finer 等价分类
%
%   out = classify_tilings_v2(A, M, reps)
%   out = classify_tilings_v2(A, M, reps, opts)
%
% 输入
%   A    : 密铺集（整数向量；内部自动取 A mod M）
%   M    : 模数
%   reps : cell 数组，reps{r} 为一个平移等价类的代表（行向量，元素在 Z_M 中，
%          满足 A + reps{r} = Z_M 为精确密铺）
%
% 两种等价关系（见论文 "Prime power tilings and tile partitions"）
%   (a) 仿射等价：B ~ u*B + t，其中 u 与 M 互素（u ∈ Z_M^×），t ∈ Z_M。
%       由 Tijdeman 膨胀引理（论文 Lemma 2.2）：A ⊕ B = Z_M 且 gcd(u,|B|)=1
%       ⇒ A ⊕ (uB) = Z_M，故仿射映像仍是 A 的补集。
%   (b) 砖块分拆等价（equivalent as tile partitions）：B ~ B' 当且仅当
%       无序砖块集合 { A+b : b ∈ B } = { A+b' : b' ∈ B' }（论文 §Partitions
%       into translates 中的定义 P_M(B) = P_M(B')）。
%
% 输出 out（结构体）
%   .M .n .Amod          模数、|A|、A mod M（升序）
%   .affine.count        仿射等价类数
%   .affine.reps         cell，每类一个代表（取该类所含平移代表中字典序最小者）
%   .affine.idx          1×numel(reps)，第 r 个平移代表所属仿射类编号
%   .affine.canon        cell，每类的典范形 min_{u,t} (uB+t)（一定是某铺法）
%   .partition.count     砖块分拆等价类数（对全部铺法逐个计算，非仅代表）
%   .partition.reps      cell，每类一个代表（该分拆类中字典序最小的铺法）
%   .partition.sizes     每个分拆类所含铺法数
%   .partitionSkipped    true/false（铺法总数过大时跳过分拆分类）
%   .skipped             '' 或 'affine'（规模超限被跳过）
%
% 选项 opts
%   .MaxOps   仿射分类的运算量上限（约 nReps*M*phi(M)*n）(默认 5e7)

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'MaxOps'), opts.MaxOps = 5e7; end

Amod = sort(mod(double(A(:))', M));
n    = numel(Amod);
if numel(unique(Amod)) ~= n
    error('classify_tilings_v2:collapse', ...
        'A mod %d 有重合元素（%d 个中只剩 %d 个），Z_%d 中不可能有密铺。', ...
        M, n, numel(unique(Amod)), M);
end

nR   = numel(reps);
R    = zeros(nR, M / n);                 % 每行一个代表（升序）
for r = 1:nR
    R(r, :) = sort(mod(double(reps{r}(:))', M));
end

% ------------------------------------------------ (a) 仿射等价
units = zeros(1, 0);
for u = 1:M - 1
    if gcd(u, M) == 1
        units(end + 1) = u;              %#ok<AGROW>
    end
end
out.skipped = '';
opsEst = double(nR) * M * numel(units) * n;   % 每个仿射像需 O(n log n) 排序
if nR > 0 && opsEst > opts.MaxOps
    out.skipped = 'affine';
    warning('classify_tilings_v2:tooMany', ...
        ['仿射分类运算量约 %g 超过 MaxOps = %g，跳过仿射分类。' ...
         '可增大 opts.MaxOps。'], opsEst, opts.MaxOps);
    out.affine = struct('count', NaN, 'reps', {{}}, 'idx', [], 'canon', {{}});
else
    key   = cell(1, nR);                 % 每个代表的仿射典范形（字符串键）
    canon = cell(1, nR);
    for r = 1:nR
        [canon{r}, key{r}] = affine_canon(R(r, :), M, units); %#ok<AGROW>
    end
    [uKeys, ~, ic] = unique(key, 'stable');
    nAff = numel(uKeys);
    aReps  = cell(1, nAff);
    aCanon = cell(1, nAff);
    aIdx   = zeros(1, nR);
    for c = 1:nAff
        members = find(ic == c);
        aIdx(members) = c;
        rows = R(members, :);
        [~, best] = sortrows(rows);      % 类内字典序最小的实际铺法
        aReps{c}  = shift_to_zero(rows(best(1), :), M);
        aCanon{c} = canon{members(1)};
    end
    out.affine = struct('count', nAff, 'reps', {sort_reps_cell(aReps)}, ...
                        'idx', aIdx, 'canon', {aCanon});
end

% ------------------------------------------------ (b) 砖块分拆等价
% 注意：分拆等价 P_M(B) = P_M(B') 不把整体平移视为等价（平移会移动全部
% 砖块），故必须对"全部铺法"（每个平移代表的所有 M 个平移）计算分拆键。
nT = M * nR;                             % 铺法总数 = numTilings（逐层已验证）
if nT * (M / n) > opts.MaxOps
    warning('classify_tilings_v2:tooMany', ...
        '分拆分类的铺法总数 %d 过大，跳过砖块分拆分类。', nT);
    out.partition = struct('count', NaN, 'reps', {{}}, 'sizes', []);
    out.partitionSkipped = true;
else
    out.partitionSkipped = false;
    pk    = cell(1, nT);
    allB  = zeros(nT, M / n);
    pos = 0;
    for r = 1:nR
        for tt = 0:M - 1
            pos = pos + 1;
            allB(pos, :) = sort(mod(R(r, :) + tt, M));
        end
    end
    for q = 1:nT
        B = allB(q, :);
        T = zeros(numel(B), n);
        for j = 1:numel(B)
            T(j, :) = sort(mod(Amod + B(j), M));
        end
        T = sortrows(T);                 % 砖块的无序集合 → 行排序后可比
        pk{q} = sprintf('%d ', T(:)');   % 整体矩阵的字符串键
    end
    % 去重放在类内做：pk 与 allB 必须保持同行对齐，不能先全局 unique
    [uPk, ~, ic] = unique(pk, 'stable');
    nPart = numel(uPk);
    pReps  = cell(1, nPart);
    pSizes = zeros(1, nPart);
    for c = 1:nPart
        members = find(ic == c);
        rows = unique(allB(members, :), 'rows');   % 去掉膨胀铺法的重复平移
        pSizes(c) = size(rows, 1);                 % 该分拆类的实际铺法数
        [~, best] = sortrows(rows);
        pReps{c} = shift_to_zero(rows(best(1), :), M);  % 类内字典序最小且 0∈B 的实际铺法
    end
    out.partition = struct('count', nPart, 'reps', {sort_reps_cell(pReps)}, ...
                           'sizes', pSizes);
end

out.M = M;  out.n = n;  out.Amod = Amod;
end

% =========================================================================
function [bestRow, bestKey] = affine_canon(B, M, units)
% B 的仿射轨道 { u*B + t } 中字典序最小的（排序后）元组
n = numel(B);
best = [];
for iu = 1:numel(units)
    Bu = sort(mod(units(iu) * B, M));    % 膨胀：gcd(u,|B|)=1 ⇒ uB 仍是补集
    for t = 0:M - 1
        row = sort(mod(Bu + t, M));
        if isempty(best) || lex_less(row, best)
            best = row;
        end
    end
end
bestRow = best;
bestKey = sprintf('%d ', best);
end

% =========================================================================
function tf = lex_less(a, b)
% 两个等长升序行向量的字典序比较（不用 sortrows，逐元素标量比较）
tf = false;
for i = 1:numel(a)
    if a(i) ~= b(i)
        tf = a(i) < b(i);
        return
    end
end
end

% =========================================================================
function B0 = shift_to_zero(B, M)
% 平移规范化：使 0 ∈ B。B-min(B) 与 B 同属一个平移类，仍是合法铺法。
B0 = sort(mod(double(B(:))' - min(double(B(:))), M));
end

% =========================================================================
function out = sort_reps_cell(reps)
% 代表列表按字典序排列（行向量等长，sortrows 逐列比较）
if isempty(reps), out = reps; return; end
nR  = numel(reps);
mat = zeros(nR, numel(reps{1}));
for r = 1:nR
    mat(r, :) = reps{r};
end
mat = sortrows(mat);
out = cell(1, nR);
for r = 1:nR
    out{r} = mat(r, :);
end
end
