%本程序为计算对数螺线的值而编写
% 绘制对数螺旋线，并添加间隔为 pi/6 的极射线
clear; clc;

% ===== 参数设置 =====
a = 2;                  % 初始半径（单位：mm）
b = 0.19;                % 螺旋参数，b = cot(螺角)
theta_max = 6 * pi;     % theta 范围
theta_step = pi/180;    % 绘图精度
delta_theta = pi/6;     % 极射线间隔
% ===== 锥角和放大系数 =====
Phi = 2*atan(b*(exp(b*2*pi)-1)/(sqrt(b^2+1)*(exp(b*2*pi)+1)));
beta = exp(b * delta_theta);
fprintf('锥角：%f', rad2deg(Phi));
fprintf('放大系数：%f', beta);
% ===== 生成螺旋数据 =====
theta = 0:theta_step:theta_max;
theta_c = 0:theta_step:theta_max-2*pi;
theta_r = 0:delta_theta:theta_max;
r = a * exp(b * theta);
r_c = 0.5 * a *(exp(b * 2 * pi) + 1) * exp(b * theta_c);
r_r = 0.5 * a *(exp(b * 2 * pi) + 1) * exp(b * theta_r);
% 极坐标转笛卡尔
x = r .* cos(theta);
y = r .* sin(theta);
x_c = r_c .* cos(theta_c);
y_c = r_c .* sin(theta_c);
x_r = r_r .* cos(theta_r);
y_r = r_r .* sin(theta_r);
% ===== 绘图设置 =====
figure;
hold on;
plot(x, y, 'b', 'LineWidth', 2);   % 对数螺旋线
plot(x_c, y_c, 'g', 'LineWidth', 2);   % 对数螺旋线中心线
scatter(x_r, y_r, 50, 'b', 'filled');  % 50个蓝色填充点，大小为50

title(['对数螺旋线，a = ', num2str(a), ', b = ', num2str(b)]);
xlabel('x (mm)');
ylabel('y (mm)');
axis equal;
grid on;

% ===== 添加极射线 =====
ray_thetas = 0:delta_theta:theta_max;
r_max = max(r);  % 射线长度设为螺旋线最大半径

for k = 1:length(ray_thetas)
    theta_k = ray_thetas(k);
    x_end = r_max * cos(theta_k);
    y_end = r_max * sin(theta_k);
    plot([0, x_end], [0, y_end], 'r--');  % 红色虚线
end

legend('对数螺旋线', '极射线');

% ===== 添加网格四边形 =====
for k = 1:length(theta_r)-1
    % 当前角度和下一个角度
    theta1 = theta_r(k);
    theta2 = theta_r(k+1);
    
    % 四个点：外螺旋上的 P1、P2，中心螺旋上的 P3、P4
    r1 = a * exp(b * theta1);
    r2 = a * exp(b * theta2);
    rc1 = 0.5 * a * (exp(b * 2*pi) + 1) * exp(b * theta1);
    rc2 = 0.5 * a * (exp(b * 2*pi) + 1) * exp(b * theta2);
    
    % 笛卡尔坐标
    x1 = r1 * cos(theta1);  y1 = r1 * sin(theta1);
    x2 = r2 * cos(theta2);  y2 = r2 * sin(theta2);
    x3 = rc2 * cos(theta2); y3 = rc2 * sin(theta2);
    x4 = rc1 * cos(theta1); y4 = rc1 * sin(theta1);
    % 大三角形三个点坐标
    A = [0, 0];
    B = [x4, y4];    %本角度的点
    C = [x3, y3];    %下一个角度的点
    
    % 边长
    AB = norm(B - A);
    BC = norm(C - B);
    CA = norm(A - C);
 
    % 角度（单位：度）
    angle_A = acosd(dot(B - A, C - A) / (norm(B - A) * norm(C - A)));
    angle_B = acosd(dot(A - B, C - B) / (norm(A - B) * norm(C - B)));
    angle_C = 180 - angle_A - angle_B;  % 或者直接用向量点乘法算
    
    % 输出
    fprintf('大三角边长:\nAB = %.3f\nBC = %.3f\nCA = %.3f\n', AB, BC, CA);
    fprintf('大三角内角 (度):\n∠A = %.2f°\n∠B = %.2f°\n∠C = %.2f°\n', ...
            angle_A, angle_B, angle_C);
     % 小三角形三个点坐标
    A = [0, 0];
    E = [x1, y1];    %本角度的点
    F = [x2, y2];    %下一个角度的点
    
    % 边长
    AE = norm(E - A);
    EF = norm(F - E);
    AF = norm(A - F);
 
    % 角度（单位：度）
    angle_A = acosd(dot(E - A, F - A) / (norm(E - A) * norm(F - A)));
    angle_E = acosd(dot(A - E, F - E) / (norm(A - E) * norm(F - E)));
    angle_F = 180 - angle_A - angle_E;  % 或者直接用向量点乘法算
    
    % 输出
    fprintf('小三角边长:\nAE = %.3f\nEF = %.3f\nAF = %.3f\n', AE, EF, AF);
    fprintf('小三角内角 (度):\n∠A = %.2f°\n∠E = %.2f°\n∠F = %.2f°\n', ...
            angle_A, angle_E, angle_F);
    % 画四边形
    fill([x1, x2, x3, x4], [y1, y2, y3, y4], [0.9 0.9 0.9], ...
     'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 1, 'HandleVisibility','off');
    % === 只处理第一个单元的几何参数 ===

end

