function PlotTimeLaw(fileName)
% Plot a trajectory. The trajectory is formatted as follows:
% time dis vel acc jerk snap
if nargin == 0
    fileName = 'traj.txt';
end
traj = load(fileName);
order = size(traj, 2) - 1;
figure('Name', 'TimeLaw')
for ii = 1:order
    subplot(order, 1, ii)
    plot(traj(:, 1), traj(:, ii+1));
end
end