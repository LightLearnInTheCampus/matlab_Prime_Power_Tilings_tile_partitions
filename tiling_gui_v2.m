function tiling_gui_v2()
%TILING_GUI_V2  单素数幂密铺枚举器的图形界面（26-11-17 版本，无绘图）
%
%   >> tiling_gui_v2
%
% 在界面中
%   1) 输入 A（元素个数为 p^alpha 的非负整数集，如 [0 3 10 13]）；
%   2) 模数 M 可留空：程序自动给出全部基周期 M = P*d（d | D）并逐一枚举；
%      也可填一个具体 M（须为 P 的倍数），只枚举 Z_M；
%   3) 点击"枚举"。输出依次为：
%        算术数据与全部基周期 → 每个基周期的平移等价类代表 →
%        仿射等价类代表 → 砖块分拆（tile partitions）等价类代表。
%   本版本不绘制任何环状彩虹图。

fig = uifigure('Name', '单素数幂密铺枚举器（基周期 × 平移/仿射/分拆）', ...
               'Position', [120 120 1000 680]);

g = uigridlayout(fig, [5 4]);
g.RowHeight    = {34, 34, 34, 34, '1x'};
g.ColumnWidth  = {96, '1x', 116, '1x'};
g.Padding      = [12 12 12 12];
g.RowSpacing   = 8;

% ---------------- 第 1 行：A ----------------
l1 = uilabel(g, 'Text', 'A =', 'FontWeight', 'bold');
l1.Layout.Row = 1; l1.Layout.Column = 1;
efA = uieditfield(g, 'text', 'Value', '[0 3 10 13]', 'FontName', 'Consolas');
efA.Layout.Row = 1; efA.Layout.Column = [2 3];
dd = uidropdown(g, 'Items', {'示例…', 'A=[0 3 10 13]  (2^2, P=4,D=5)', ...
                             'A=[0 3 6]      (3^1, P=9)', ...
                             'A=[0 10 20]    (3^1, P=3,D=10)', ...
                             'A=[0 4]        (2^1, P=8)', ...
                             'A=[0 1 2 3 4 5 6 7]  (2^3)'}, ...
                'Value', '示例…', 'ValueChangedFcn', @(~, ~) fill_example());
dd.Layout.Row = 1; dd.Layout.Column = 4;

% ---------------- 第 2 行：M ----------------
l2 = uilabel(g, 'Text', 'M =');
l2.Layout.Row = 2; l2.Layout.Column = 1;
efM = uieditfield(g, 'text', 'Value', '', ...
                  'Placeholder', '留空 → 枚举全部基周期 P*d (d|D)', ...
                  'FontName', 'Consolas');
efM.Layout.Row = 2; efM.Layout.Column = 2;
l2b = uilabel(g, 'Text', '（填 M：须为 P 的倍数，只枚举 Z_M）', ...
              'FontColor', [0.4 0.4 0.4]);
l2b.Layout.Row = 2; l2b.Layout.Column = [3 4];

% ---------------- 第 3 行：选项 ----------------
l3 = uilabel(g, 'Text', '枚举上限');
l3.Layout.Row = 3; l3.Layout.Column = 1;
efMax = uieditfield(g, 'numeric', 'Value', 200000, 'Limits', [1 1e8], ...
                    'RoundFractionalValues', 'on');
efMax.Layout.Row = 3; efMax.Layout.Column = 2;
cb3 = uicheckbox(g, 'Text', '约化到基周期（A 或补集为陪集并时缩小）', 'Value', false);
cb3.Layout.Row = 3; cb3.Layout.Column = [3 4];

% ---------------- 第 4 行：按钮 ----------------
btn = uibutton(g, 'Text', '枚举', 'FontWeight', 'bold', ...
               'ButtonPushedFcn', @onRun);
btn.Layout.Row = 4; btn.Layout.Column = [1 2];
btn2 = uibutton(g, 'Text', '清空输出', 'ButtonPushedFcn', @(~, ~) clear_output());
btn2.Layout.Row = 4; btn2.Layout.Column = [3 4];

% ---------------- 第 5 行：输出 ----------------
ta = uitextarea(g, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 12);
ta.Layout.Row = 5; ta.Layout.Column = [1 4];

set(ta, 'Value', { ...
    '说明：输入 A（|A| = p^alpha）。模数 M 留空时程序给出全部基周期', ...
    '      M = P*d（d | D），对每个基周期枚举 Z_M 的全部铺法并给出：', ...
    '      平移等价类代表 → 仿射等价类代表 → 砖块分拆等价类代表。', ...
    '      M 可以比 A 中的元素小（程序自动取 A mod M）。', ...
    '      本版本不绘制环状彩虹图。'});

    % ===================== 回调 =====================
    function onRun(~, ~)
        try
            A = str2num(efA.Value); %#ok<ST2NM>
            if isempty(A)
                error('无法解析 A，请写成 [0 3 10 13] 这样的形式。');
            end
            if isempty(strtrim(efM.Value))
                M = [];
            else
                M = str2num(efM.Value); %#ok<ST2NM>
                if isempty(M) || ~isscalar(M)
                    error('无法解析 M，请填正整数或留空。');
                end
            end
            opts = struct('ReducePeriod', cb3.Value, ...
                          'MaxEnum', efMax.Value, ...
                          'ShowReps', 40);
            txt = tiling_report_v2(A, M, opts);
            set(ta, 'Value', strsplit(txt, newline));
        catch ME
            set(ta, 'Value', strsplit(sprintf('运行出错：%s\n\n%s', ME.message, ...
                getReport(ME, 'basic', 'hyperlinks', 'off')), newline));
        end
    end

    function clear_output()
        set(ta, 'Value', {''});
    end

    function fill_example()
        switch dd.Value
            case 'A=[0 3 10 13]  (2^2, P=4,D=5)'
                efA.Value = '[0 3 10 13]';  efM.Value = '';
            case 'A=[0 3 6]      (3^1, P=9)'
                efA.Value = '[0 3 6]';      efM.Value = '9';
            case 'A=[0 10 20]    (3^1, P=3,D=10)'
                efA.Value = '[0 10 20]';    efM.Value = '';
            case 'A=[0 4]        (2^1, P=8)'
                efA.Value = '[0 4]';        efM.Value = '';
            case 'A=[0 1 2 3 4 5 6 7]  (2^3)'
                efA.Value = '[0 1 2 3 4 5 6 7]'; efM.Value = '';
        end
    end
end
