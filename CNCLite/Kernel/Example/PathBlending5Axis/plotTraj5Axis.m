%% load data.
traj = load('traj5Axis.txt');
wcs = load('wcsData.txt');
mcs = load('mcsData.txt');

linearLimits = [50, 200, 2000]; % in mm.
angularLimits = [5, 50, 500]; % in degree.

%% data process
ts = traj(2, 1);
aLinear = zeros(size(traj, 1), 1);
jLinear = aLinear;
aAngular = aLinear;
jAngular = aLinear;
aLinear(2:end) = diff(traj(:, 3), 1, 1) / ts;
jLinear(2:end) = diff(aLinear, 1, 1) / ts;
aAngular(2:end) = diff(traj(:, 8), 1, 1) / ts;
jAngular(2:end) = diff(aAngular, 1, 1) / ts;

pos = traj(:, 18:22); % axial position
vel = zeros(size(pos) );
acc = vel;
vel(2:end, :) = diff(pos, 1, 1) / ts; % axial velocity
acc(2:end, :) = diff(vel, 1, 1) / ts; % axial acceleration

%% draw figures
wl = 2;
cl = {'r', 'g', 'b', 'c', 'm'};

%% Axial position
figure('Name', 'Axial position')
plot(traj(:, 1), pos(:, 1), 'r', traj(:, 1), pos(:, 2), 'g', ...
    traj(:, 1), pos(:, 3), 'b', traj(:, 1), pos(:, 4), 'c', ...
    traj(:, 1), pos(:, 5), 'm', 'LineWidth', wl);
ylabel('各轴位置 (mm(deg))');
xlabel('时间(s)');
legend('X', 'Y', 'Z', 'A', 'C');
set(gca, 'FontName', '黑体');

%% Axial velocity
figure('Name', 'Axial velocity')
plot(traj(:, 1), vel(:, 1), 'r', traj(:, 1), vel(:, 2), 'g', ...
    traj(:, 1), vel(:, 3), 'b', traj(:, 1), vel(:, 4), 'c', ...
    traj(:, 1), vel(:, 5), 'm', 'LineWidth', wl);
ylabel('各轴速度 (mm(deg)/s)');
xlabel('时间 (s)');
legend('X', 'Y', 'Z', 'A', 'C');
set(gca, 'FontName', '黑体');

%% Axial acceleration
% figure('Name', 'Axial acceleration')
% plot(traj(:, 1), acc(:, 1), 'r', traj(:, 1), acc(:, 2), 'g', ...
%     traj(:, 1), acc(:, 3), 'b', traj(:, 1), acc(:, 4), 'c', ...
%     traj(:, 1), acc(:, 5), 'm', 'LineWidth', wl);
% ylabel('各轴加速度 (mm(deg)/s^2)');
% xlabel('时间 (s)');
% legend('X', 'Y', 'Z', 'A', 'C');
% set(gca, 'FontName', '黑体');

%% % tool positions in WCS
interval = 1; % plot the tool path every interval points.
toolLength = 20;
colorLow = [0.0, 0.0, 1.0];
colorHigh = [1.0, 0.0, 0.0];
feedLow = min(traj(:, 3) );
feedHigh = max(traj(:, 3) );
feedDelta = feedHigh - feedLow;
figure('Name', 'Position Contour');
plot3(traj(:, 12), traj(:, 13), traj(:, 14) ); % position curve.
hold on;
for i = 1 : size(wcs, 1)
    plot3([wcs(i, 1), wcs(i, 1) + wcs(i, 4)*toolLength],...
        [wcs(i, 2), wcs(i, 2) + wcs(i, 5)*toolLength],...
        [wcs(i, 3), wcs(i, 3) + wcs(i, 6)*toolLength],...
        'k-o', 'LineWidth', wl, 'MarkerSize', 2);
end
numColor = 256;
% color map
cmap = [linspace(colorLow(1), colorHigh(1), numColor)', linspace(colorLow(2), colorHigh(2), numColor)',...
    linspace(colorLow(3), colorHigh(3), numColor)'];
for i = 1 : size(traj, 1)
    if mod(i, interval) == 0
        alpha = (traj(i, 3) - feedLow) / feedDelta;
        indexColor = round(alpha*numColor) + 1;
        if indexColor > numColor
            indexColor = numColor;
        end
        col = cmap(indexColor, :);
        plot3([traj(i, 12), traj(i, 12) + traj(i, 15) * toolLength], ...
            [traj(i, 13), traj(i, 13) + traj(i, 16) * toolLength], ...
            [traj(i, 14), traj(i, 14) + traj(i, 17) * toolLength], 'Color', col);
        hold on;
    end
end
hold off;
axis equal;
shading flat;
cb = colorbar('Limits', [feedLow, feedHigh]);
cb.Label.String = '{\bff}(mm/s)';
cb.Label.FontName = '黑体';
cb.Label.FontSize = 10;
caxis([feedLow, feedHigh]); % colormap map the value between [feedLow, feedHigh] linearly by colormap
set(gcf, 'Colormap',cmap);
xlabel('{\bfx}(mm)');
ylabel('{\bfy}(mm)');
zlabel('{\bfz}(mm)');
set(gca, 'FontName', '黑体');

%% tool orientation
interval = 1; % plot the tool path every interval points.
toolLength = 20;
colorLow = [0.0, 0.0, 1.0];
colorHigh = [1.0, 0.0, 0.0];
feedLow = min(traj(:, 8) );
feedHigh = max(traj(:, 8) );
feedDelta = feedHigh - feedLow;

% plot as 2-D figure
figure('Name', 'Orientation Contour');
% tool orientation
plot(wcs(:, 4), wcs(:, 5), 'ko', 'MarkerSize', 4 );
hold on;
numColor = 256;
% color map
cmap = [linspace(colorLow(1), colorHigh(1), numColor)', linspace(colorLow(2), colorHigh(2), numColor)',...
    linspace(colorLow(3), colorHigh(3), numColor)'];
for i = 1 : size(traj, 1)-1
    alpha = (traj(i, 8) - feedLow) / feedDelta;
    indexColor = round(alpha*numColor) + 1;
    if indexColor > numColor
        indexColor = numColor;
    end
    col = cmap(indexColor, :);
    plot([traj(i, 15), traj(i+1, 15)], ...
        [traj(i, 16), traj(i+1, 16)],... 
        'Color', col, 'LineWidth', wl);
    hold on;
end
% plot the iso-curve for k-component.
num_iso = 200;
k_iso = 0.9704 : 0.006 : 0.9944;
% k_iso = 0.92 : 0.01 : 0.97;
for ii = 1 : length(k_iso)
    i_max = sqrt( 1 - k_iso(ii)^2 );
    i_iso = linspace(-i_max, i_max, num_iso);
    j_iso = sqrt(i_max^2 - i_iso.^2);
    plot(i_iso, j_iso, 'k-.');
    plot(i_iso, -j_iso, 'k-.');
end
hold off;
axis equal;
xlim([-0.2, 0.05]);
ylim([-0.22, 0.1]);
% xlim([-0.36, -0.07]);
% ylim([-0.2, 0]);
% xlim([-0.36, -0.24]);
% ylim([-0.2, -0.11]);
cb = colorbar('Limits', [feedLow, feedHigh]);
cb.Label.String = '{\omega}(deg/s)';
cb.Label.FontName = '黑体';
cb.Label.FontSize = 10;
caxis([feedLow, feedHigh]); % colormap map the value between [feedLow, feedHigh] linearly by colormap
set(gcf, 'Colormap',cmap);
xlabel('{\bfi}');
ylabel('{\bfj}');
zlabel('{\bfk}');
set(gca, 'FontName', '黑体');




%% linear velocity
figure('Name', 'Linear velocities');
plot(traj(:, 1), traj(:, 3), 'r-', 'LineWidth', wl);
hold on;
% upper and lower bounds
plot([traj(1, 1), traj(end, 1)], linearLimits(1)*[1, 1], 'k-.');
hold off;
ylabel('刀尖点速度 (mm/s)');
xlabel('时间 (s)');
set(gca, 'FontName', '黑体');

%% angular velocity
figure('Name', 'Angular velocities');
plot(traj(:, 1), traj(:, 8), 'r-', 'LineWidth', wl);
hold on;
% upper and lower bounds
plot([traj(1, 1), traj(end, 1)], angularLimits(1)*[1, 1], 'k-.');
hold off;
ylabel('刀轴转速 (deg/s)');
xlabel('时间 (s)');
set(gca, 'FontName', '黑体');

%% linear acceleration
figure('Name', 'Linear acceleration');
plot(traj(:, 1), aLinear, 'r-', 'LineWidth', wl);
hold on;
% upper and lower bounds
plot([traj(1, 1), traj(end, 1)], linearLimits(2)*[1, 1], 'k-.');
plot([traj(1, 1), traj(end, 1)], -linearLimits(2)*[1, 1], 'k-.');
hold off;
ylabel('刀尖点加速度(mm/s^2)');
xlabel('时间 (s)');
set(gca, 'FontName', '黑体');

%% linear jerk
figure('Name', 'Linear jerk');
plot(traj(:, 1), jLinear, 'r-', 'LineWidth', wl);
hold on;
% upper and lower bounds
plot([traj(1, 1), traj(end, 1)], linearLimits(3)*[1, 1], 'k-.');
plot([traj(1, 1), traj(end, 1)], -linearLimits(3)*[1, 1], 'k-.');
hold off;
ylabel('刀尖点跃度 (mm/s^3)');
xlabel('时间 (s)');
set(gca, 'FontName', '黑体');

%% angular acceleration
figure('Name', 'Angular acceleration');
plot(traj(:, 1), aAngular, 'r-', 'LineWidth', wl);
hold on;
% upper and lower bounds
plot([traj(1, 1), traj(end, 1)], angularLimits(2)*[1, 1], 'k-.');
plot([traj(1, 1), traj(end, 1)], -angularLimits(2)*[1, 1], 'k-.');
hold off;
ylabel('刀轴转动加速度(deg/s^2)');
xlabel('时间 (s)');
set(gca, 'FontName', '黑体');

%% angular jerk
figure('Name', 'Angular jerk');
plot(traj(:, 1), jAngular, 'r-', 'LineWidth', wl);
hold on;
% upper and lower bounds
plot([traj(1, 1), traj(end, 1)], angularLimits(3)*[1, 1], 'k-.');
plot([traj(1, 1), traj(end, 1)], -angularLimits(3)*[1, 1], 'k-.');
hold off;
ylabel('刀轴转动跃度 (deg/s^3)');
xlabel('时间 (s)');
set(gca, 'FontName', '黑体');
