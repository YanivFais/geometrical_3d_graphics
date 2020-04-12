% Weighted iterative closest point algorithm based on kd tree fast nearest
% neighbor query.
%
% (C) Alex Bronstein, 2008. All rights reserved.
function [surface_y, R, t, dist_y, nangle_y, dist_x, nangle_x] = icp(surface_y, surface_x, weight_y, weight_x, options)

% Settings
if nargin < 5,
    options = [];
end

if nargin == 3,
   options = weight_y; 
end

if nargin < 4,
   weight_x = [];
   weight_y = [];    
end

if ~isfield(options,'verbose'),                 verbose = 1;                 else    verbose = options.verbose;                           end
if ~isfield(options,'min_normal_angle'),        min_normal_angle = 0.5;      else    min_normal_angle = options.min_normal_angle;         end
if ~isfield(options,'max_vertex_dist'),         
    % Automatically estimate the rejection threshold as 10% of the shape
    % size
    X = [surface_x.X surface_x.Y surface_x.Z];
    x0 = mean(X, 1);
    X0 = X - repmat(x0, [size(X,1) 1]);
    max_vertex_dist  = 0.1*sqrt(sum(sum(X0.^2))/3);     
else
    max_vertex_dist  = options.max_vertex_dist;    
end
if ~isfield(options,'max_iter'),                max_iter = 1000;             else    max_iter = options.max_iter;                         end
if ~isfield(options,'min_coord_change'),        min_coord_change = 1e-3;     else    min_coord_change = options.min_coord_change;         end
if ~isfield(options,'min_rel_error_change'),    min_rel_error_change = 1e-5; else    min_rel_error_change = options.min_rel_error_change; end
if ~isfield(options,'initial_centering'),       initial_centering = 1;       else    initial_centering = 0;                               end

if isempty(weight_x),
    weight_x = ones(length(surface_x.X), 1);
end
if isempty(weight_y),
    weight_y = ones(length(surface_y.X), 1);
end

% Initialization
TRANS = eye(4,4);

if verbose,
    fprintf(1, ' Starting iterative closest point algorithm.\n');
    fprintf(1, ' Iter \t        Error \t  Rel. change \t  Coord change \t Points \n');
end

% Compute surface normals
[Nv,Nt] = compute_normals (surface_x);
surface_x.Nv = Nv;
surface_x.Nt = Nt;

[Nv,Nt] = compute_normals (surface_y);
surface_y.Nv = Nv;
surface_y.Nt = Nt;

% Bring centroid into correspondence
if initial_centering,
    X = [surface_x.X surface_x.Y surface_x.Z];
    x0 = sum(X .* repmat(weight_x(:)/sum(weight_x), [1 3]), 1);
    Y = [surface_y.X surface_y.Y surface_y.Z];
    y0 = sum(Y .* repmat(weight_y(:)/sum(weight_y), [1 3]),1);
    Y = Y' + repmat(x0(:)-y0(:), [1 size(Y,1)]);
    surface_y.X = Y(1,:)';
    surface_y.Y = Y(2,:)';
    surface_y.Z = Y(3,:)';
end

% Build nearest neighbor kd-tree
if ~isfield(options,'tree_x'),
    tree_x = ann('init', [surface_x.X surface_x.Y surface_x.Z]');
else
    tree_x = options.tree_x;
end

error_history = [];

for iter = 1:max_iter,
     
    if verbose > 1,
        figure(1);
        subplot(1,2,1);
        hx = trisurf(surface_x.TRIV,surface_x.X,surface_x.Y,surface_x.Z); axis image; axis off;
        hold on;
        hy = trisurf(surface_y.TRIV,surface_y.X,surface_y.Y,surface_y.Z);
        hold off;
        shading interp;
        set(hx,'FaceColor',[1 0 0]);
        set(hy,'FaceColor',[0 0 1]);
        lighting phong;
        camlight head;
        title(sprintf('Iteration %d', iter));
        drawnow;
    end
    
    % Find closest points
    [corr_x, dist] = ann('search', tree_x, [surface_y.X surface_y.Y surface_y.Z]', 1, 'eps', 1);
    
    % Compute angles between corresponding normals
    nx = surface_x.Nv(corr_x,:);
    nangle = sum(nx.*surface_y.Nv, 2);

    % Reject bad correspondences
    corr_y = find(nangle(:) >= min_normal_angle & dist(:) <= max_vertex_dist);

    % If correspondence is empty, take 50% of the closest points
    if isempty(corr_y),
       [dist, corr_y] = sort(dist);
       corr_y = corr_y(1:round(length(corr_y)*0.5));
    end
    
    corr_x = corr_x(corr_y);
    

    % Compute weights
    w = weight_x(corr_x) .* weight_y(corr_y);
    w = w/sum(w);

    % Solve for best rigid transformation using SVD    
    X = [surface_x.X(corr_x) surface_x.Y(corr_x) surface_x.Z(corr_x)];
    x0 = sum(X .* repmat(w(:), [1 3]), 1);
    X0 = X - repmat(x0, [size(X,1) 1]);

    Y = [surface_y.X(corr_y) surface_y.Y(corr_y) surface_y.Z(corr_y)];
    y0 = sum(Y .* repmat(w(:), [1 3]),1);
    Y0 = Y - repmat(y0, [size(Y,1) 1]);

    H = Y0'*spdiags(w(:), 0, length(w), length(w))*X0;
    [U,S,V] = svd(H');
    R = U*V';
    t = x0(:) - R*y0(:);
    
    % Compute error for first iteration
    if iter == 1,
        Y = [surface_y.X(corr_y) surface_y.Y(corr_y) surface_y.Z(corr_y)];
        X = [surface_x.X(corr_x) surface_x.Y(corr_x) surface_x.Z(corr_x)];
        err = sqrt(sum(sum((X-Y).^2.*repmat(w(:), [1 3]))));
    end
    
    % Update the target surface
    Y = [surface_y.X surface_y.Y surface_y.Z];
    Y_old = Y';
    Y = R*Y' + repmat(t(:), [1 size(Y,1)]);
    surface_y.X = Y(1,:)';
    surface_y.Y = Y(2,:)';
    surface_y.Z = Y(3,:)';
    surface_y.Nv = Nv*R';
    surface_y.Nt = Nt*R';
    
    % Cumulative transformation
    TRANS = [R t; 0 0 0 1]* TRANS;
    
    % Change in coordinates
    coord_change = sqrt(sum(sum((Y-Y_old).^2)));
    
    % Compute error
    Y = [surface_y.X(corr_y) surface_y.Y(corr_y) surface_y.Z(corr_y)];
    X = [surface_x.X(corr_x) surface_x.Y(corr_x) surface_x.Z(corr_x)];
    err_old = err;
    err = sqrt(sum(sum((X-Y).^2.*repmat(w(:), [1 3]))));
    rel_err_change = (err_old-err)/err_old;
    
    error_history(iter) = err;
    
    if verbose > 1,
        figure(1);
        subplot(1,2,2);
        plot(error_history);
        xlabel('Iteration');
        ylabel('Alignment error');
        drawnow;
    end    
        
    % Print iteration information
    if verbose,
        fprintf(1,' %4d \t %12.6g \t %12.2g \t  %12.6g \t %4.1f%%\n', ...
                iter, err, rel_err_change, coord_change, length(corr_x)/length(surface_x.X)*100);
    end
    
    if rel_err_change >= 0 & abs(rel_err_change) < min_rel_error_change,
        if verbose,
            fprintf(1, ' Minimum reached: relative error change %g < %g.\n\n', rel_err_change, min_rel_error_change);
        end
        break;
    end

    if coord_change < min_coord_change,
        if verbose,
            fprintf(1, ' Minimum reached: coordinate change %g < %g.\n\n', coord_change, min_coord_change);
        end
        break;
    end
        
end

if iter == max_iter,
    if verbose,
        fprintf(1, ' Maximum number of %d iterations reached. Optimization stopped.\n\n', max_iter);
    end
end


% Compute distances and normal angles
if nargout >= 4,
    [corr_x, dist_y] = ann('search', tree_x, [surface_y.X surface_y.Y surface_y.Z]', 1, 'eps', 1);
    nx = surface_x.Nv(corr_x,:);
    nangle_y = sum(nx.*surface_y.Nv, 2);
    dist_y = dist_y(:);
end

if nargout >= 6,
    if ~isfield(options,'tree_y'),
        tree_y = ann('init', [surface_y.X surface_y.Y surface_y.Z]');
    else
        tree_y = options.tree_y;
    end
    [corr_y, dist_x] = ann('search', tree_y, [surface_x.X surface_x.Y surface_x.Z]', 1, 'eps', 1);
    ny = surface_y.Nv(corr_y,:);
    nangle_x = sum(ny.*surface_x.Nv, 2);
    if ~isfield(options,'tree_y'),
        ann('deinit', tree_y);
    end
    dist_x = dist_x(:);
end

% Extract transformation
R = TRANS(1:3,1:3);
t = TRANS(1:3,4);

% Done using the kd tree
if ~isfield(options,'tree_x'),
    ann('deinit', tree_x);
end

