%DEMO_TILING_V2  26-11-17 版本命令行演示：全部基周期 + 三种等价类代表
%   运行：>> demo_tiling_v2
%   等价于在 GUI（tiling_gui_v2）中输入相应 A 后点击"枚举"。

fprintf('########## 例 1：A = [0 3 10 13]（论文算例，全部基周期） ##########\n\n');
txt = tiling_report_v2([0 3 10 13]);
fprintf('%s\n\n', txt);

fprintf('########## 例 2：A = [0 1 18 19]（基周期 {4, 12, 36}） ##########\n\n');
txt = tiling_report_v2([0 1 18 19], [], struct('ShowReps', 12));
fprintf('%s\n\n', txt);

fprintf('########## 例 3：A = [0 3 10 13]，只看 Z_20，约化到基周期 ##########\n\n');
txt = tiling_report_v2([0 3 10 13], 20, struct('ReducePeriod', true, 'ShowReps', 10));
fprintf('%s\n', txt);
