% Numerical geometry of nonrigid shapes.
% (C) Alexander & Michael Bronstein, 2008.
% All rights reserved.

close all;
figure(1);
%set(gcf, 'Position',[1 36 1600 1085]);
pos = get(gcf, 'Position');
pos([1 3]) = [1 1600];
set(gcf, 'Position',pos);
clf;
set(gcf, 'Color', [1 1 1]);

rand('seed', 0);

% Load shapes
load michael3 shape; shape_x = shape;
load michael4 shape; shape_y = shape;

% Run ICP - successful convergence
opt = set_options('max_vertex_dist', Inf, 'min_normal_angle', -1, ...
    'max_iter', 250, 'initial_centering', 1, 'verbose', 2, ...
    'min_rel_error_change', 1e-6, 'min_coord_change', 1e-4);
wx = ones(size(shape_x.X));
wy = ones(size(shape_y.X));
[shape_y_transformed, R, t] = icp(shape_y, shape_x, wy, wx, opt);

fprintf(1, 'Press any key to continue...\n');
pause;

figure(1);
set(gcf, 'Position',pos);
clf;
set(gcf, 'Color', [1 1 1]);

rand('seed', 0);

% Load shapes
load victoria1 shape; shape_x = shape;

% Run ICP - successful convergence
opt = set_options('max_vertex_dist', Inf, 'min_normal_angle', -1, ...
    'max_iter', 250, 'initial_centering', 1, 'verbose', 2, ...
    'min_rel_error_change', 1e-6, 'min_coord_change', 1e-4);
wx = ones(size(shape_x.X));
wy = ones(size(shape_y.X));
[shape_y_transformed, R, t] = icp(shape_y, shape_x, wy, wx, opt);

