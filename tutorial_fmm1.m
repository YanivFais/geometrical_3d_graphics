% Numerical geometry of nonrigid shapes.
% (C) Alexander & Michael Bronstein, 2008.
% All rights reserved.

rand('seed', 0);
load michael3;

for k=1:10,
    src = round(rand*(length(shape.X)-1))+1;
    D0 = repmat(Inf, [length(shape.X) 1]);
    D0(src) = 0;
    D = fastmarch(shape.TRIV, shape.X, shape.Y, shape.Z, D0, struct('mode', 'single'));

    
    trisurf(shape.TRIV, shape.X, shape.Y, shape.Z, D); 
    axis image; axis off; shading interp; lighting phong; 
    view([-10 15]); camlight head; colormap hot;
    title(sprintf('Source at vertex %d', src));
    drawnow;
    
    pause;
end
