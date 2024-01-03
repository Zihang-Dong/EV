% 初始化YALMIP和Gurobi
clear;
clc;
close all;
yalmip('clear');

N = 4; % 电动汽车数量
T = 24; % 时间窗口数量（小时）

% 定义参数
dt = 1; % 采样时间间隔 (小时)
P_max = 8; % 充放电功率上界 (kW)
E_max = 40; % 电池最大容量 (kWh)
E_min = 0; % 电池最小容量 (kWh)
E_initial = 20*rand(N, 1); % 初始电量 (kWh)
E_final_min = 35; % 最终最小能量要求 (kWh)
eta_c = 0.95; % 充电效率
eta_d = 0.92; % 放电效率

% 分时电价设置（7小时低价，其余高价）
low_price = 0.385; % 波谷电价 (￥/kWh)
high_price = 0.725; % 波峰电价 (￥/kWh)
buy_price = high_price * ones(T, 1);
buy_price(1:7) = low_price; % 低价时段为1-7点和22-24点
buy_price(22:24) = low_price;
sell_price = 0.75*buy_price; % 卖电价格 = 买电价格*0.6

% 定义优化变量
P_c = sdpvar(T, N); % 充电功率
P_d = sdpvar(T, N); % 放电功率
z_c = binvar(T, N); % 充电状态二元变量
z_d = binvar(T, N); % 放电状态二元变量
E = sdpvar(T+1, N); % 电池能量

% 目标函数：最小化电费成本
Objective = 0;
for i = 1:N
    for t = 1:T
        Objective = Objective + buy_price(t) * P_c(t, i) + sell_price(t) * P_d(t, i);
    end
end

% 约束条件
Constraints = [];
for i = 1:N
    % 赋值初始电量
    E(1, i) = E_initial(i);

    for t = 1:T
        % 功率约束
        Constraints = [Constraints, -P_max*z_d(t, i) <= P_d(t, i), P_d(t, i) <= 0];
        Constraints = [Constraints, 0 <= P_c(t, i), P_c(t, i) <= P_max*z_c(t, i)];
        
        % 电池能量动态变化
        Constraints = [Constraints, E(t+1, i) == 0.999*E(t, i) + (eta_c*P_c(t, i) + 1/eta_d*P_d(t, i))*dt];
        
        % 电池能量边界
        Constraints = [Constraints, E_min <= E(t, i), E(t, i) <= E_max];
        
        % 状态互斥约束（一辆车同时只能充电或放电）
        Constraints = [Constraints, z_c(t, i) + z_d(t, i) <= 1];
        
        % 电池最终最小能量要求
        if t == T
            Constraints = [Constraints, E(t+1, i) >= E_final_min];
        end
    end
end

% 求解优化问题
options = sdpsettings('verbose', 1, 'solver', 'gurobi');
sol = optimize(Constraints, Objective, options);

% 检查结果
if sol.problem == 0
    % 成功求解
    P_optimal = value(P_c+P_d);
    E_optimal = value(E);
    Z_c_optimal = value(z_c);
    Z_d_optimal = value(z_d);
    disp('Optimal charging/discharging schedule:');
    disp(P_optimal);
else
    disp('Something went wrong!');
    disp(sol.info);
end


%% 图

% 绘制电价
figure;
plot(0:T-1, buy_price, '-o', 'LineWidth', 1.5);
hold on
plot(0:T-1, sell_price, '-*', 'LineWidth', 1.5);
title('Electricity Price over 24 Hours');
xlabel('Time (hours)');
ylabel('Price (￥/kWh)');
grid on;
xlim([0, T]);
legend('Buy', 'Sell', 'Location', 'best');

% 进行绘制
EVs_to_plot = 1:N;

% 绘制充放电功率曲线和电池能量变化曲线
figure;
for i = 1:length(EVs_to_plot)
    EV_idx = EVs_to_plot(i);

    % 绘制充放电功率
    subplot(4, 2, 2*i-1);
    plot(0:T-1, P_optimal(:, EV_idx), 'LineWidth', 1.5);
    title(['EV ' num2str(EV_idx) ' Charging/Discharging Power']);
    xlabel('Time (hours)');
    ylabel('Power (kW)');
    grid on;
    xlim([0, T]);
    ylim([-P_max, P_max]);
    yticks([-8 -4 0 4 8]);

    % 绘制电池能量变化曲线
    subplot(4, 2, 2*i);
    plot(0:T, E_optimal(:, EV_idx), 'LineWidth', 1.5);
    title(['EV ' num2str(EV_idx) ' Battery Energy']);
    xlabel('Time (hours)');
    ylabel('Energy (kWh)');
    grid on;
    xlim([0, T]);
    ylim([0, E_max]);
    yticks([0 10 20 30 40]);
end

