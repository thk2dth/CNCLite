traj = load('trajS2.txt'); % load data
ts = traj(2, 1); % sampling period.
pos = traj(:, 7:9);
vel = zeros(size(pos) );
vel(2:end, :) = diff(pos, 1, 1) / ts;

%% Draw figures
wl = 2;
%% Trajectory
figure('Name', 'Trajectory');
subplot(4, 1, 1)
plot(traj(:, 1), traj(:, 3), 'LineWidth', wl);
ylabel('{\bfV} (deg/s)');
subplot(4, 1, 2)
plot(traj(:, 1), traj(:, 4), 'LineWidth', wl);
ylabel('{\bfA} (deg/s^2)');
subplot(4, 1, 3)
plot(traj(:, 1), traj(:, 5), 'LineWidth', wl);
ylabel('{\bfJ} (deg/s^3)');
subplot(4, 1, 4)
plot(traj(:, 1), traj(:, 6), 'LineWidth', wl);
ylabel('{\bfS} (deg/s^4)');
xlabel('{\bfT} (s)');
set(gca, 'FontName', 'Times New Roman');

%% Contour
figure('Name', 'Contour');
hold on;
sphere(256);
shading interp;
alpha(0.5);
plot3(pos(:, 1), pos(:, 2), pos(:, 3), 'k', 'LineWidth', wl);
xlabel('{\bfI}');
ylabel('{\bfJ}');
zlabel('{\bfK}');
hold off;
axis equal;
set(gca, 'FontName', 'Times New Roman');

%% Axial velocity
figure('Name', 'Axial velocity');
plot(traj(:, 1), vel(:, 1), 'r', traj(:, 1), vel(:, 2), 'g');
hold on;
plot(traj(:, 1), traj(:, 3), 'k', 'LineWidth', wl);
legend('X', 'Y', 'T');
xlabel('{\bfT} (s)');
ylabel('{\bfV} (deg/s)');
set(gca, 'FontName', 'Times New Roman');