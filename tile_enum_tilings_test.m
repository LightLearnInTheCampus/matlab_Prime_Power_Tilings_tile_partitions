function tile_enum_tilings_test()
%TILE_ALL_TILINGS_TEST  用论文算例与穷举法校验 tile_enum_tilings
%
% 校验内容
%  (1) 论文 §7 的全部数值结果（P, D, F_A, E_A, 平移类数）
%  (2) 与"穷举所有子集 + 唯一性检验"的暴力算法逐一对拍（小规模）
%  (3) 每个输出代表确实满足 A ⊕ B = Z_M
%
% 直接运行： >> tile_enum_tilings_test

clc
fprintf('=============== 论文算例校验 ===============\n');
fprintf('实际加载：%s\n', which('tile_enum_tilings'));

% ---------- 例 1：A = {0,3,10,13} ----------
A1 = [0 3 10 13];
o  = tile_enum_tilings(A1, [], struct('Report', false));
fprintf('文件版本：%s\n', o.version);
assert(o.isTile);
assert(isequal(o.nu, [1 2]) && o.p == 2 && o.k == 2);
assert(isequal(o.Delta, [1 5]) && o.D == 5 && o.P == 4 && o.PD == 20);
assert(isequal(o.canonicalTile, [0 1 10 11]));
assert(o.table(1,3) == 4  && o.table(1,4) == 4);      % F_A(4)=4,  E_A(4)=4
assert(o.table(2,3) == 64 && o.table(2,4) == 60);     % F_A(20)=64, E_A(20)=60
assert(o.numClasses == 4);                             % 1 + 3
fprintf('例 1  A=[0 3 10 13] :  P=%d D=%d, F(4)=%d E(4)=%d, F(20)=%d E(20)=%d, 类数=%d  [OK]\n', ...
    o.P, o.D, o.table(1,3), o.table(1,4), o.table(2,3), o.table(2,4), o.numClasses);

% 论文列出的三个基本周期 20 的代表（应与枚举结果同属三类）
ref = {[0 2 4 6 8], [0 4 6 8 12], [0 2 6 8 14]};
for r = 1:3
    assert(any_in_same_class(ref{r}, o.representatives, 20), ...
        '论文代表 B%d = %s 不出现在枚举结果的任何平移类中；程序代表为 %s', ...
        r, mat2str(ref{r}), rep_str(o.representatives));
end
fprintf('       论文 B1,B2,B3 均出现在枚举结果中  [OK]\n');

% ---------- 例 2：A = {0,3,6} ----------
A2 = [0 3 6];
o2 = tile_enum_tilings(A2, 9, struct('Report', false));
assert(o2.P == 9 && o2.D == 1 && o2.numTilings == 27 && o2.numClasses == 3);
ref2 = {[0 1 2], [0 1 5], [0 2 4]};
for r = 1:3
    assert(any_in_same_class(ref2{r}, o2.representatives, 9), ...
        '论文代表 %s 不出现在枚举结果的任何平移类中；程序代表为 %s', ...
        mat2str(ref2{r}), rep_str(o2.representatives));
end
fprintf('例 2  A=[0 3 6]     :  P=%d, 补集=%d, 平移类=%d, 三个代表吻合  [OK]\n', ...
    o2.P, o2.numTilings, o2.numClasses);

% ---------- 例 3：A = {0,10,20} ----------
A3 = [0 10 20];
o3 = tile_enum_tilings(A3, [], struct('Report', false));
assert(o3.P == 3 && o3.D == 10 && o3.PD == 30);
expF = [3 9 243 59049];
expE = [3 6 240 58800];
expC = [1 1 16 1960];
assert(isequal(o3.table(:,3)', expF), 'F_A 不符');
assert(isequal(o3.table(:,4)', expE), 'E_A 不符');
assert(isequal(o3.table(:,5)', expC), '平移类数不符');
assert(o3.numClasses == 1978);
fprintf('例 3  A=[0 10 20]   :  F=%s  E=%s  类=%d（论文 1978）  [OK]\n', ...
    mat2str(expF), mat2str(expE), o3.numClasses);

% ---------- 例 4：A = {0,4}（nu_k > k 的情形）----------
o4 = tile_enum_tilings([0 4], [], struct('Report', false));
assert(o4.nu == 3 && o4.P == 8 && o4.D == 1);
assert(o4.numTilings == 16 && o4.numClasses == 2);
fprintf('例 4  A=[0 4]       :  nu=%d, P=%d, 补集=%d, 类=%d  [OK]\n', ...
    o4.nu, o4.P, o4.numTilings, o4.numClasses);

% ---------- 例 5：A = {0,1,...,7} ----------
o5 = tile_enum_tilings(0:7, [], struct('Report', false));
assert(isequal(o5.nu, [1 2 3]) && o5.P == 8 && o5.numTilings == 8 && o5.numClasses == 1);
fprintf('例 5  A=[0..7]      :  nu=[1 2 3], 补集=%d, 类=%d  [OK]\n', ...
    o5.numTilings, o5.numClasses);

% ---------- 例 6：非密铺集 ----------
o6 = tile_enum_tilings([0 1 3], [], struct('Report', false));
assert(~o6.isTile);
fprintf('例 6  A=[0 1 3]     :  正确判为"非密铺集"  [OK]\n');

% ---------- 例 7：不可达模数 ----------
o7 = tile_enum_tilings([0 3 6], 10, struct('Report', false));
assert(o7.numTilings == 0 && isempty(o7.representatives));
fprintf('例 7  A=[0 3 6], M=10 : P=9 不整除 10，正确报告"无铺法"  [OK]\n');

% ================== 与暴力算法对拍 ==================
fprintf('\n============ 与暴力穷举对拍 ============\n');
cases = { [0 1],      3;        % (A, M)
          [0 1],      4;
          [0 4],      8;
          [0 3],      4;
          [0 1 2],    3;
          [0 1 2],    6;
          [0 3 6],    9;
          [0 1 2 3],  4;
          [0 1 2 3],  8;
          [0 2 4 6],  8;
          [0 3 10 13], 20;
          [0 1 4 5],  8;
          [0 2 5 7],  16 };
allok = true;
for c = 1:size(cases,1)
    A = cases{c,1};  M = cases{c,2};
    o = tile_enum_tilings(A, M, struct('Report', false, 'MaxEnum', 1e6));
    [nBf, kBF] = brute_force(A, M);
    kOur = cellfun(@(B) canon_key(B, M), o.representatives, 'UniformOutput', false);
    ok1 = (nBf == o.numTilings);
    ok2 = (numel(kBF) == o.numClasses) && all(ismember(kBF, kOur)) && all(ismember(kOur, kBF));
    % 每个代表必须真正密铺
    ok3 = true;
    for r = 1:numel(o.representatives)
        if ~is_tiling(A, o.representatives{r}, M), ok3 = false; end
    end
    allok = allok && ok1 && ok2 && ok3;
    fprintf('A=%-14s M=%-3d  铺法: 本程序 %-6d 暴力 %-6d %s | 类: %-4d %-4d %s | 代表有效 %s\n', ...
        mat2str(A), M, o.numTilings, nBf, tick(ok1), o.numClasses, numel(kBF), ...
        tick(ok2), tick(ok3));
end
% ---------------- 基周期约化（A 或补集为子群陪集之并）的校验 ----------------
fprintf('\n============ 基周期约化校验 ============\n');
rcases = { [0 3 6], 9 ; [0 4], 8 ; [0 1], 8 ; [0 1 2 3], 8 ; [0 3 10 13], 20 };
rok = true;
for c = 1:size(rcases, 1)
    A = rcases{c, 1};  M = rcases{c, 2};
    A0 = A - min(A);
    [A2, M2, ri] = reduce_tiling_period(A0, M);
    okA = true;  note = '';
    if M2 < M
        o2  = tile_enum_tilings(A2, M2, struct('Report', false));
        okA = o2.isTile && o2.numClasses >= 1;            % 约化后仍能密铺 Z_M*
        % 计数关系：只发生 A 型约化时 F_A(Z_M) = m^{|B|} * F_A(Z_M*)
        if okA && ~any(strcmp(ri.kinds, 'B'))
            o1  = tile_enum_tilings(A, M, struct('Report', false));
            okA = (o1.numTilings == ri.nTilingsFactor * o2.numTilings);
            note = sprintf('（F: %d = %d × %d）', o1.numTilings, ...
                ri.nTilingsFactor, o2.numTilings);
        end
        % 可靠性：约化后的代表按"B 型步取完全提升"回抬，应给出 Z_M 的铺法
        if okA
            B = o2.representatives{1};
            for st = numel(ri.kinds):-1:1
                if strcmp(ri.kinds{st}, 'B')
                    d = ri.dmods(st + 1);
                    B = lift_set(B, d, ri.dmods(st) / d);
                end
            end
            okA = is_tiling(A0, mod(B, M), M);
            if ~okA, note = '（回抬后不是 Z_M 的铺法）'; end
        end
    else
        note = '（无需约化）';
    end
    rok = rok && okA;
    fprintf('A=%-14s M=%-3d → M*=%-3d 倍数 %-2d %s %s\n', ...
        mat2str(A), M, M2, ri.factor, tick(okA), note);
end
if ~rok, error('基周期约化校验失败。'); end
fprintf('约化校验通过。\n');

if allok
    fprintf('\n全部对拍通过。\n');
else
    error('对拍失败。');
end
end

% -------------------------------------------------------------------------
function C = lift_set(B, d, m)
% 把 B ⊂ Z_d 完全提升为 {b + k*d : k = 0..m-1} ⊂ Z_{d*m}
C = [];
for k = 0:m-1
    C = [C, B(:)' + k * d]; %#ok<AGROW>
end
C = mod(unique(C), d * m);
end

% -------------------------------------------------------------------------
function s = tick(ok)if ok, s = '[OK]'; else, s = '[X ]'; end
end

% -------------------------------------------------------------------------
function tf = is_tiling(A, B, M)
% A ⊕ B = Z_M 的精确检验（用累加计数，避免 A(:)+B(:).' 的隐式扩展）
A = mod(A(:)', M); B = mod(B(:)', M);
if numel(A) * numel(B) ~= M, tf = false; return; end
cnt = zeros(M, 1);
for i = 1:numel(A)
    cnt = cnt + accumarray(mod(A(i) + B, M)' + 1, 1, [M 1]);
end
tf = all(cnt == 1);
end

% -------------------------------------------------------------------------
function [n, keys] = brute_force(A, M)
% 穷举所有 |B| = M/|A| 的子集 B ⊂ Z_M，收集真正密铺的平移类
nA = numel(A);
if mod(M, nA) ~= 0, n = 0; keys = {}; return; end
t = M / nA;
if t < 1 || t > M, n = 0; keys = {}; return; end
C = nchoosek(0:M-1, t);
n = 0; keys = {};
for r = 1:size(C, 1)
    if is_tiling(A, C(r, :), M)
        n = n + 1;
        k = canon_key(C(r, :), M);
        if ~any(strcmp(k, keys)), keys{end+1} = k; end %#ok<AGROW>
    end
end
end

% -------------------------------------------------------------------------
function k = canon_key(B, M)
% 平移等价的不变量：轨道上所有集合"排序后元组"的字典序最小者（用 sortrows，不做字符比较）
B = sort(unique(mod(B(:)', M)));
t = numel(B);
if t == 0, k = ''; return; end
R = zeros(M, t);                       % M 个平移，每行是该平移的排序元组
for s = 1:M
    R(s, :) = sort(mod(B + s - 1, M));
end
R = sortrows(R, 1:t);                  % 按第 1 列到第 t 列依次比较
k = mat2str(R(1, :));
end

% -------------------------------------------------------------------------
function tf = same_class(B, C, M)
% B 与 C 是否循环平移等价：直接比对轨道，不借助任何"键"
B = sort(mod(B(:)', M));
C = sort(mod(C(:)', M));
tf = false;
if numel(B) ~= numel(C), return; end
for s = 0:M-1
    if isequal(sort(mod(B + s, M)), C), tf = true; return; end
end
end

% -------------------------------------------------------------------------
function tf = any_in_same_class(B, reps, M)
% B 的平移类是否等于 reps 中某一项所在的类
tf = false;
for j = 1:numel(reps)
    if same_class(B, reps{j}, M), tf = true; return; end
end
end

% -------------------------------------------------------------------------
function s = rep_str(reps)
% 把代表列表拼成一行，便于断言失败时直接看出程序给了什么
s = '';
for j = 1:numel(reps)
    s = [s mat2str(reps{j}) '  ']; %#ok<AGROW>
end
if isempty(strtrim(s)), s = '(空)'; end
s = strtrim(s);
end

