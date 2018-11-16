traj = load('traj3Axis.txt'); % load data
ts = traj(2, 1); % sampling period.
pos = traj(:, 7:9);
vel = zeros(size(pos) );
vel(2:end, :) = diff(pos, 1, 1) / ts;

% original point
pts = load('../../../NC Program/OriginalSpiralDataTest.dat');

%% Draw figures
wl = 2;
figure('Name', 'Trajectory');
subplot(4, 1, 1)
plot(traj(:, 1), traj(:, 3), 'LineWidth', wl);
ylabel('{\bfV} (mm/s)');
subplot(4, 1, 2)
plot(traj(:, 1), traj(:, 4), 'LineWidth', wl);
ylabel('{\bfA} (mm/s^2)');
subplot(4, 1, 3)
plot(traj(:, 1), traj(:, 5), 'LineWidth', wl);
ylabel('{\bfJ} (mm/s^3)');
subplot(4, 1, 4)
plot(traj(:, 1), traj(:, 6), 'LineWidth', wl);
ylabel('{\bfS} (mm/s^4)');
xlabel('{\bfT} (s)');
set(gca, 'FontName', 'Times New Roman');

figure('Name', 'Contour');
plot(pos(:, 1), pos(:, 2), 'b', 'LineWidth', wl);
hold on;
plot(pts(:, 1), pts(:, 2), 'k');
xlabel('{\bfX} (mm)');
ylabel('{\bfY} (mm)');
set(gca, 'FontName', 'Times New Roman');

figure('Name', 'Axial velocity');
plot(traj(:, 1), vel(:, 1), 'r', traj(:, 1), vel(:, 2), 'g');
hold on;
plot(traj(:, 1), traj(:, 3), 'k', 'LineWidth', wl);
legend('X', 'Y', 'T');
xlabel('{\bfT} (s)');
ylabel('{\bfV} (mm/s)');
set(gca, 'FontName', 'Times New Roman');

figure('Name', 'Curvature')
plot(traj(:, 1), traj(:, 11) );
xlabel('{\bf u}');
ylabel('{\bf \kappa}{\it(mm^{-1})}');