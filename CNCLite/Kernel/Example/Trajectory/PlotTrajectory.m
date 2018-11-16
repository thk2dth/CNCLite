function PlotTrajectory(fileName)
% Plot a trajectory. The trajectory is formatted as follows:
% time dis vel acc jerk snap
if nargin == 0
    fileName = 'traj.txt';
end
traj = load(fileName);
% the last six columns are axial positions and axial velocities.
order = size(traj, 2) - 1 - 6; 
dim = 3;
figure('Name', 'TimeLaw')
for ii = 1:order
    subplot(order, 1, ii)
    plot(traj(:, 1), traj(:, ii+1));
end

figure('Name', 'Contour')
plot3(traj(:, order+2), traj(:, order+3), traj(:, order+4));

figure('Name', 'AxialVelocity')
for ii = 1:dim
    subplot(dim, 1, ii)
    plot(traj(:, 1), traj(:, order+1+dim+ii));
end

end