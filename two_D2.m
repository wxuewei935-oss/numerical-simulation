% 清除环境
clear; close all; clc;

%% 捕食-被捕食模型参数
a = 1.0;        % 猎物的自然增长率
b = 0.5;        % 捕食导致的猎物死亡率  
c = 1.0;       % 捕食者转化效率
d = 1.0;        % 捕食者死亡率
e = 1.0;
L = 2;         % 空间域长度 [0, L]
N = floor(L*pi/0.2);        % 空间网格数 (将有 N+1 个点)
T = 500;        % 模拟总时间

D_u = 0.8;        % u物种的扩散系数
D_v = 1.0;        % v物种的扩散系数

% 计算空间步长和网格点
dx = L*pi / N;
x = linspace(0, L*pi, N+1); % 空间网格

% 时间跨度
tspan = [0 T];

%% 初始条件 (u和v物种空间分离)
% u物种: 
u_0 = 2*ones(N+1,1)'+0.001*rand(N+1,1)';
% v物种: 
v_0 = 2*ones(N+1,1)'+0.001*rand(N+1,1)';

u0_system = [u_0; v_0];

%% 求解捕食-被捕食系统
[t_data, u_data] = ode15s(@(t, u) predatorPreySystem(t, u, D_u, D_v, a, b, c, d, e, dx, N), tspan, u0_system);

u_species = u_data(:, 1:N+1);      % u物种数据
v_species = u_data(:, N+2:end);    % v物种数据

%% 全局变量声明（用于回调函数）
global h_bar_u h_bar_v h_final_bar_u h_final_bar_v ...
       h_time_line1 h_time_line2 ...
       u_species_data v_species_data t_species_data x_data L_data ...
       time_slider time_value progress_label ...
       play_button play_timer is_playing speed_slider bar_width ...
       plot_threshold threshold_slider threshold_value  % 添加绘图阈值相关变量

% 存储原始数据
u_species_data = u_species;
v_species_data = v_species;
t_species_data = t_data;
x_data = x;
L_data = L;

% 设置绘图阈值 - 小于此值的在柱状图中绘制为0
plot_threshold = 0.001;  % 默认阈值

% 计算最大值（使用原始数据）
max_u = max(u_species(:));
max_v = max(v_species(:));
max_combined = max(max_u, max_v) * 1.1;

% 初始时间索引
initial_time_idx = 1;

% 柱状图宽度调整
bar_width = 0.35 * (x(2)-x(1)); % 减少宽度，留出间隙

%% 辅助函数：应用阈值到数据
function data_processed = applyThreshold(data, threshold)
    % 将小于阈值的数据设置为0
    data_processed = data;
    data_processed(abs(data_processed) < threshold) = 0;
end

%% 创建主图窗口（调整布局确保交互控件不挡住图形）
fig = figure('Position', [50, 50, 1600, 950], 'Name', 'u-v物种系统可视化 (带阈值处理)', ...
    'NumberTitle', 'off', 'Resize', 'on');

% 设置图形窗口的布局，为底部控件留出空间
% 使用 normalized 单位确保布局适配窗口大小

%% 子图1: u物种密度分布热图
ax1 = subplot(3, 4, [1, 2, 5, 6]);
set(ax1, 'Position', [0.05, 0.40, 0.44, 0.55]); % 调整位置
h_heat_u = pcolor(x, t_data, u_species);
shading interp; 
colormap(ax1, jet);
colorbar('Position', [0.495, 0.40, 0.02, 0.55]);
xlabel('Spatial position x', 'FontSize', 11);
ylabel('time', 'FontSize', 11);
title('Thermogram of u species density distribution', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold on;
% 添加当前时间线
h_time_line1 = plot([min(x), max(x)], [t_data(initial_time_idx), t_data(initial_time_idx)], 'w-', 'LineWidth', 2.5);

%% 子图2: v物种密度分布热图
ax2 = subplot(3, 4, [3, 4, 7, 8]);
set(ax2, 'Position', [0.55, 0.40, 0.44, 0.55]); % 调整位置
h_heat_v = pcolor(x, t_data, v_species);
shading interp; 
colormap(ax2, hot);
colorbar('Position', [1.00, 0.40, 0.02, 0.55]);
xlabel('Spatial position x', 'FontSize', 11);
ylabel('time', 'FontSize', 11);
title('Thermogram of v species density distribution', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold on;
% 添加当前时间线
h_time_line2 = plot([min(x), max(x)], [t_data(initial_time_idx), t_data(initial_time_idx)], 'w-', 'LineWidth', 2.5);

%% 子图3: 当前时刻分布柱状图（并排显示）
global ax3;
ax3 = subplot(3, 4, [9, 10]);
set(ax3, 'Position', [0.05, 0.15, 0.44, 0.22]); % 调整位置

% 获取初始时刻数据并应用阈值
initial_u_raw = u_species(initial_time_idx, :);
initial_v_raw = v_species(initial_time_idx, :);
initial_u_processed = applyThreshold(initial_u_raw, plot_threshold);
initial_v_processed = applyThreshold(initial_v_raw, plot_threshold);

% 统计处理后数据中非零值的数量
u_nonzero = sum(initial_u_processed ~= 0);
v_nonzero = sum(initial_v_processed ~= 0);

% 设置图例文本
u_legend_text = sprintf('u (non zero point: %d/%d)', u_nonzero, length(initial_u_raw));
v_legend_text = sprintf('v (non zero point: %d/%d)', v_nonzero, length(initial_v_raw));

% 计算偏移量，使两个柱状图并排显示
bar_offset = bar_width * 0.5;

% 创建并排的柱状图
% u物种柱状图（左侧）
bar_positions_u = x - bar_offset;
h_bar_u = bar(bar_positions_u, initial_u_processed, bar_width, ...
    'FaceColor', [0.2, 0.8, 0.2], 'EdgeColor', [0, 0.6, 0], ...
    'LineWidth', 1, 'DisplayName', u_legend_text);

hold on;

% v物种柱状图（右侧）
bar_positions_v = x + bar_offset;
h_bar_v = bar(bar_positions_v, initial_v_processed, bar_width, ...
    'FaceColor', [0.9, 0.2, 0.2], 'EdgeColor', [0.7, 0, 0], ...
    'LineWidth', 1, 'DisplayName', v_legend_text);

hold off;
xlabel('Spatial position x', 'FontSize', 11);
ylabel('Density', 'FontSize', 11);
title(sprintf('Current time distribution - histogram (t = %.1f, threshold=%.4f)', t_data(initial_time_idx), plot_threshold), ...
    'FontSize', 12, 'FontWeight', 'bold');
legend({u_legend_text, v_legend_text}, 'FontSize', 9, 'Location', 'northeast');
grid on;

% 动态调整y轴范围
current_max = max(max(initial_u_processed), max(initial_v_processed));
if current_max == 0
    current_max = 0.1; % 如果全部为0，设置一个默认范围
end
axis([min(x)-dx max(x)+dx 0 current_max*1.3]);

%% 子图4: 最终时刻分布柱状图
ax4 = subplot(3, 4, [11, 12]);
set(ax4, 'Position', [0.55, 0.15, 0.44, 0.22]); % 调整位置
final_time_idx = length(t_data);

% 获取最终时刻数据并应用阈值
final_u_raw = u_species(final_time_idx, :);
final_v_raw = v_species(final_time_idx, :);
final_u_processed = applyThreshold(final_u_raw, plot_threshold);
final_v_processed = applyThreshold(final_v_raw, plot_threshold);

% 统计处理后数据中非零值的数量
u_final_nonzero = sum(final_u_processed ~= 0);
v_final_nonzero = sum(final_v_processed ~= 0);

% 设置图例文本
u_final_legend = sprintf('u (non zero point: %d/%d)', u_final_nonzero, length(final_u_raw));
v_final_legend = sprintf('v (non zero point: %d/%d)', v_final_nonzero, length(final_v_raw));

% 创建并排的柱状图
% 最终时刻u物种柱状图（左侧）
h_final_bar_u = bar(bar_positions_u, final_u_processed, bar_width, ...
    'FaceColor', [0.2, 0.8, 0.2], 'EdgeColor', [0, 0.6, 0], ...
    'LineWidth', 1, 'DisplayName', u_final_legend);

hold on;

% 最终时刻v物种柱状图（右侧）
h_final_bar_v = bar(bar_positions_v, final_v_processed, bar_width, ...
    'FaceColor', [0.9, 0.2, 0.2], 'EdgeColor', [0.7, 0, 0], ...
    'LineWidth', 1, 'DisplayName', v_final_legend);

hold off;
xlabel('Spatial position x', 'FontSize', 11);
ylabel('Density', 'FontSize', 11);
title(sprintf('Final time distribution - histogram (t = %.1f, threshold=%.4f)', t_data(final_time_idx), plot_threshold), ...
    'FontSize', 12, 'FontWeight', 'bold');
legend({u_final_legend, v_final_legend}, 'FontSize', 9, 'Location', 'northeast');
grid on;

% 动态调整y轴范围
final_max = max(max(final_u_processed), max(final_v_processed));
if final_max == 0
    final_max = 0.1; % 如果全部为0，设置一个默认范围
end
axis([min(x)-dx max(x)+dx 0 final_max*1.3]);

%% 添加控制按钮面板（放在更底部）
control_panel = uipanel('Position', [0.05, 0.02, 0.90, 0.05], ...
    'BackgroundColor', [0.9 0.9 0.9], ...
    'BorderType', 'line', 'HighlightColor', [0.5 0.5 0.5]);

% 创建播放/暂停按钮
play_button = uicontrol('Parent', control_panel, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.02, 0.1, 0.08, 0.8], ...
    'String', '▶ 播放', ...
    'FontSize', 9, ...
    'Callback', @togglePlay);

% 创建步进按钮
step_back_button = uicontrol('Parent', control_panel, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.11, 0.1, 0.07, 0.8], ...
    'String', '⏮ 后退', ...
    'FontSize', 9, ...
    'Callback', {@stepTime, -1});

step_forward_button = uicontrol('Parent', control_panel, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.19, 0.1, 0.07, 0.8], ...
    'String', '⏭ 前进', ...
    'FontSize', 9, ...
    'Callback', {@stepTime, 1});

% 创建重置按钮
reset_button = uicontrol('Parent', control_panel, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.27, 0.1, 0.07, 0.8], ...
    'String', '↺ 重置', ...
    'FontSize', 9, ...
    'Callback', @resetTime);

% 创建速度控制标签
speed_label = uicontrol('Parent', control_panel, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.35, 0.1, 0.05, 0.8], ...
    'String', '速度:', ...
    'FontSize', 9, ...
    'HorizontalAlignment', 'right', ...
    'BackgroundColor', [0.9 0.9 0.9]);

speed_slider = uicontrol('Parent', control_panel, ...
    'Style', 'slider', ...
    'Units', 'normalized', ...
    'Position', [0.41, 0.1, 0.12, 0.8], ...
    'Min', 1, ...
    'Max', 10, ...
    'Value', 5, ...
    'SliderStep', [0.1, 0.2]);

export_gif_button = uicontrol('Parent', control_panel, ...
    'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [0.82, 0.1, 0.12, 0.8], ...  % 位置根据实际调整
    'String', '导出柱状图GIF', ...
    'FontSize', 9, ...
    'Callback', @exportBarGIF);

% 创建绘图阈值控制标签
threshold_label = uicontrol('Parent', control_panel, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.54, 0.1, 0.05, 0.8], ...
    'String', '绘图阈值:', ...
    'FontSize', 9, ...
    'HorizontalAlignment', 'right', ...
    'BackgroundColor', [0.9 0.9 0.9]);

threshold_slider = uicontrol('Parent', control_panel, ...
    'Style', 'slider', ...
    'Units', 'normalized', ...
    'Position', [0.60, 0.1, 0.12, 0.8], ...
    'Min', 0.0001, ...
    'Max', 0.1, ...
    'Value', plot_threshold, ...
    'SliderStep', [0.01, 0.05], ...
    'Callback', @updateThreshold);

% 创建阈值显示标签
threshold_value = uicontrol('Parent', control_panel, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.73, 0.1, 0.08, 0.8], ...
    'String', sprintf('%.4f', plot_threshold), ...
    'FontSize', 9, ...
    'FontWeight', 'bold', ...
    'BackgroundColor', [0.9 0.9 0.9]);

%% 添加时间轴滑块面板（放在底部）
slider_panel = uipanel('Position', [0.05, 0.00, 0.90, 0.02], ...
    'BackgroundColor', [0.85 0.85 0.85], ...
    'BorderType', 'none');

% 创建时间滑块
time_slider = uicontrol('Parent', slider_panel, ...
    'Style', 'slider', ...
    'Units', 'normalized', ...
    'Position', [0.0, 0.0, 1.0, 1.0], ...
    'Min', 1, ...
    'Max', length(t_data), ...
    'Value', initial_time_idx, ...
    'SliderStep', [1/(length(t_data)-1), 10/(length(t_data)-1)], ...
    'Callback', @updateTime);

% 创建时间显示和进度标签（放在滑块上方）
time_display_panel = uipanel('Position', [0.05, 0.07, 0.90, 0.02], ...
    'BackgroundColor', [0.92 0.92 0.92], ...
    'BorderType', 'none');

% 时间显示标签
time_label = uicontrol('Parent', time_display_panel, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.0, 0.0, 0.1, 1.0], ...
    'String', '当前时间:', ...
    'FontSize', 9, ...
    'BackgroundColor', [0.92 0.92 0.92]);

time_value = uicontrol('Parent', time_display_panel, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.1, 0.0, 0.1, 1.0], ...
    'String', sprintf('%.1f', t_data(initial_time_idx)), ...
    'FontSize', 9, ...
    'FontWeight', 'bold', ...
    'BackgroundColor', [0.92 0.92 0.92]);

progress_label = uicontrol('Parent', time_display_panel, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.8, 0.0, 0.2, 1.0], ...
    'String', sprintf('进度: %.1f%%', 100*initial_time_idx/length(t_data)), ...
    'FontSize', 9, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', ...
    'BackgroundColor', [0.92 0.92 0.92]);

%% 初始化播放状态变量
is_playing = false;
play_timer = timer('ExecutionMode', 'fixedRate', 'Period', 0.1, ...
    'TimerFcn', @playTimerCallback, 'BusyMode', 'queue');

%% 回调函数定义

% 更新时间显示
function updateTime(source, ~)
    global u_species_data v_species_data t_species_data x_data L_data ...
           h_bar_u h_bar_v h_final_bar_u h_final_bar_v ...
           h_time_line1 h_time_line2 time_value progress_label ...
           bar_width plot_threshold  % 添加绘图阈值
    
    time_idx = round(source.Value);
    t_current = t_species_data(time_idx);
    
    % 更新滑块值显示
    time_value.String = sprintf('%.1f', t_current);
    
    % 更新当前时刻分布柱状图
    axes_handle = get(h_bar_u, 'Parent'); % 获取axes句柄
    
    % 获取当前时刻原始数据并应用阈值
    current_u_raw = u_species_data(time_idx, :);
    current_v_raw = v_species_data(time_idx, :);
    
    % 应用阈值处理
    current_u_processed = applyThreshold(current_u_raw, plot_threshold);
    current_v_processed = applyThreshold(current_v_raw, plot_threshold);
    
    % 统计非零值数量
    u_nonzero = sum(current_u_processed ~= 0);
    v_nonzero = sum(current_v_processed ~= 0);
    
    % 计算偏移位置
    bar_offset = bar_width * 0.5;
    bar_positions_u = x_data - bar_offset;
    bar_positions_v = x_data + bar_offset;
    
    % 更新柱状图数据（使用处理后的数据）
    set(h_bar_u, 'XData', bar_positions_u, 'YData', current_u_processed);
    set(h_bar_v, 'XData', bar_positions_v, 'YData', current_v_processed);
    
    % 更新子图标题（显示阈值信息）
    title(axes_handle, sprintf('Current time distribution - histogram (t = %.1f, threshold=%.4f))', t_current, plot_threshold), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    % 更新图例（显示非零点数量）
    u_legend_text = sprintf('u (non zero point: %d/%d)', u_nonzero, length(current_u_raw));
    v_legend_text = sprintf('v (non zero point: %d/%d)', v_nonzero, length(current_v_raw));
    legend(axes_handle, {u_legend_text, v_legend_text}, 'FontSize', 9, 'Location', 'northeast');
    
    % 动态调整y轴范围
    current_max = max(max(current_u_processed), max(current_v_processed));
    if current_max == 0
        current_max = 0.1; % 如果全部为0，设置一个默认范围
    end
    dx_val = x_data(2) - x_data(1);
    axis(axes_handle, [min(x_data)-dx_val max(x_data)+dx_val 0 current_max*1.3]);
    
    % 更新热图中的时间线
    set(h_time_line1, 'YData', [t_current, t_current]);
    set(h_time_line2, 'YData', [t_current, t_current]);
    
    % 更新进度信息
    progress_label.String = sprintf('进度: %.1f%%', 100*time_idx/length(t_species_data));
    
    % 刷新图形
    drawnow;
end

% 更新阈值函数
function updateThreshold(source, ~)
    global plot_threshold threshold_value time_slider h_final_bar_u h_final_bar_v ...
           u_species_data v_species_data t_species_data x_data bar_width
    
    % 更新阈值
    plot_threshold = source.Value;
    
    % 更新阈值显示
    threshold_value.String = sprintf('%.4f', plot_threshold);
    
    % 更新最终时刻柱状图（重新应用阈值）
    final_time_idx = length(t_species_data);
    final_u_raw = u_species_data(final_time_idx, :);
    final_v_raw = v_species_data(final_time_idx, :);
    
    % 应用阈值处理
    final_u_processed = applyThreshold(final_u_raw, plot_threshold);
    final_v_processed = applyThreshold(final_v_raw, plot_threshold);
    
    % 统计非零值数量
    u_final_nonzero = sum(final_u_processed ~= 0);
    v_final_nonzero = sum(final_v_processed ~= 0);
    
    % 更新最终时刻柱状图数据
    bar_offset = bar_width * 0.5;
    bar_positions_u = x_data - bar_offset;
    bar_positions_v = x_data + bar_offset;
    
    set(h_final_bar_u, 'XData', bar_positions_u, 'YData', final_u_processed);
    set(h_final_bar_v, 'XData', bar_positions_v, 'YData', final_v_processed);
    
    % 更新最终时刻图例
    u_final_legend = sprintf('u (non zero point: %d/%d)', u_final_nonzero, length(final_u_raw));
    v_final_legend = sprintf('v (non zero point: %d/%d)', v_final_nonzero, length(final_v_raw));
    
    final_axes = get(h_final_bar_u, 'Parent');
    legend(final_axes, {u_final_legend, v_final_legend}, 'FontSize', 9, 'Location', 'northeast');
    
    % 更新最终时刻标题
    title(final_axes, sprintf('最终时刻分布 - 柱状图 (t = %.1f, 阈值=%.4f)', ...
        t_species_data(final_time_idx), plot_threshold), 'FontSize', 12, 'FontWeight', 'bold');
    
    % 调整最终时刻y轴范围
    final_max = max(max(final_u_processed), max(final_v_processed));
    if final_max == 0
        final_max = 0.1; % 如果全部为0，设置一个默认范围
    end
    dx_val = x_data(2) - x_data(1);
    axis(final_axes, [min(x_data)-dx_val max(x_data)+dx_val 0 final_max*1.3]);
    
    % 更新当前时刻显示
    updateTime(time_slider, []);
end

% 播放/暂停切换
function togglePlay(~, ~)
    global is_playing play_timer play_button
    
    is_playing = ~is_playing;
    if is_playing
        play_button.String = '⏸ 暂停';
        start(play_timer);
    else
        play_button.String = '▶ 播放';
        stop(play_timer);
    end
end

% 步进函数
function stepTime(~, ~, step)
    global time_slider t_species_data
    
    new_idx = round(time_slider.Value) + step;
    new_idx = max(1, min(new_idx, length(t_species_data)));
    time_slider.Value = new_idx;
    updateTime(time_slider, []);
end

% 重置时间
function resetTime(~, ~)
    global time_slider is_playing
    
    time_slider.Value = 1;
    updateTime(time_slider, []);
    if is_playing
        togglePlay([], []);
    end
end

% 播放定时器回调
function playTimerCallback(~, ~)
    global time_slider speed_slider t_species_data is_playing
    
    speed = get(speed_slider, 'Value');
    step_size = round(speed / 5); % 根据速度调整步长
    
    new_idx = round(time_slider.Value) + step_size;
    if new_idx > length(t_species_data)
        new_idx = length(t_species_data);
        togglePlay([], []); % 到达末尾时停止
    end
    
    time_slider.Value = new_idx;
    updateTime(time_slider, []);
end

function exportBarGIF(~, ~)
    global t_species_data time_slider ax3
    global is_playing play_timer
    
    % ----- 用户可调参数 -----
    gif_filename = 'bar_animation_slow.gif';   % 输出文件名
    frame_step = 2;          % 每隔2帧保存一帧（数值越大文件越小）
    delay_time = 0.3;        % 每帧显示0.3秒（约3.3帧/秒，较慢）
    margin_pixels = 25;      % 裁剪边距（像素），确保轴标签和标题完整
    % -----------------------
    
    if isempty(ax3) || ~isvalid(ax3)
        errordlg('未找到柱状图坐标轴，请确保 ax3 已定义为全局变量', '导出失败');
        return;
    end
    
    % 暂停播放定时器
    if exist('play_timer', 'var') && isa(play_timer, 'timer') && strcmp(play_timer.Running, 'on')
        stop(play_timer);
        was_playing = true;
        is_playing = false;
    else
        was_playing = false;
    end
    
    num_steps = length(t_species_data);
    fprintf('开始生成柱状图 GIF（保留坐标轴刻度和名称），共 %d 帧（实际保存 %d 帧）...\n', ...
        num_steps, ceil(num_steps/frame_step));
    
    fig = gcf;
    set(fig, 'Color', 'w');  % 白色背景，更清晰
    
    % 获取图形窗口和子图的像素位置（用于精确裁剪）
    oldFigUnits = get(fig, 'Units');
    set(fig, 'Units', 'pixels');
    figPos = get(fig, 'Position');
    set(fig, 'Units', oldFigUnits);
    
    oldAxUnits = get(ax3, 'Units');
    set(ax3, 'Units', 'pixels');
    axPos = get(ax3, 'Position');   % [left bottom width height] 相对于 figure 左下角
    set(ax3, 'Units', oldAxUnits);
    
    % 计算裁剪区域：在子图四周增加边距，确保坐标轴标签、标题、图例完整
    crop_x = max(1, axPos(1) - margin_pixels);
    crop_y = max(1, figPos(4) - (axPos(2) + axPos(4) + margin_pixels));
    crop_w = min(figPos(3) - crop_x, axPos(3) + 2*margin_pixels);
    crop_h = min(figPos(4) - crop_y, axPos(4) + 2*margin_pixels);
    crop_rect = [crop_x-30, crop_y, crop_w+30, crop_h+15];
    
    h_waitbar = waitbar(0, '正在生成 GIF，请稍候...');
    
    for idx = 1:frame_step:num_steps
        % 更新时间滑块，触发柱状图更新
        time_slider.Value = idx;
        updateTime(time_slider, []);
        drawnow;
        
        % 捕获整个图形窗口并裁剪出子图区域（包含坐标轴装饰）
        frame = getframe(fig);
        img = frame2im(frame);
        img_cropped = imcrop(img, crop_rect);
        [imind, cm] = rgb2ind(img_cropped, 256);
        
        % 写入 GIF
        if idx == 1
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end
        
        waitbar(idx/num_steps, h_waitbar);
    end
    
    close(h_waitbar);
    fprintf('✅ GIF 已保存为：%s\n', gif_filename);
    
    % 恢复播放状态
    if was_playing
        start(play_timer);
        is_playing = true;
    end
end

%% 捕食-被捕食模型函数
function dudt = predatorPreySystem(~, u, D_u, D_v, a, b, c, d, e, dx, N)
    u_species = u(1:N+1);
    v_species = u(N+2:end);
    
    dudt_u = zeros(N+1, 1);
    dudt_v = zeros(N+1, 1);
    
    % u物种的扩散
    for i = 2:N
        dudt_u(i) = D_u * (u_species(i-1) - 2*u_species(i) + u_species(i+1)) / dx^2;
    end
    dudt_u(1) = D_u * (-1*u_species(1) + u_species(2)) / dx^2;
    dudt_u(N+1) = D_u * (u_species(N) - 1*u_species(N+1)) / dx^2;
    
    % v物种的扩散
    for i = 2:N
        dudt_v(i) = D_v * (v_species(i-1) - 2*v_species(i) + v_species(i+1)) / dx^2;
    end
    dudt_v(1) = D_v * (-2*v_species(1) + v_species(2)) / dx^2;
    dudt_v(N+1) = D_v * (v_species(N) - 1*v_species(N+1)) / dx^2;
    
    % 相互作用反应项
    % u物种: u_species.*(a+b * u_species) - c * u_species .* v_species
    % v物种: d * u_species .* v_species - e * v_species.*v_species
    dudt_u = dudt_u + u_species.*(a+b * u_species) - c * u_species .* v_species;
    dudt_v = dudt_v + d * u_species .* v_species - e * v_species.*v_species;
    
    dudt = [dudt_u; dudt_v];
end