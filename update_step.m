function [ d_3 ] = update_step( V,d_ , Velocity , normal )
%fast marching update step from two vertices and distances 
% (ignoring obtuse triangles)
% V in the form (P1,P2) where P1/2 in the form (x;y)
% d = [d1;d2]
% Veclocity is velocity input in (x1,x2,x3)
% normal the normal over the surface at x1,x2

% need to solve:
% d_3^2*(1^T*Q*1) - d_3*(2*1^T*Q*d) + d^T*Q*d-1=0

d = d_ .* [Velocity(1);Velocity(2)];

Q = V'*V;

if (d(1)==Inf)
  if (d(2)==Inf)
      d_3 = Inf;
      return;
  else 
      d_3 = normal(2)*Velocity(2)+d_(2); % distance from P2
  end
else
  if (d(2)==Inf)
      d_3 = normal(1)*Velocity(1)+d_(1); % distance from P1
  end

  a = [1,1]*Q*[1;1];
  b = -2*[1,1]*Q*d;
  c = d'*Q*d-1;
  [d_3a,d_3b ] = quadratic(a,b,c);
  d_3 = max(d_3a,d_3b);
end

if (~isreal(d_3) | d_3<0)
    d_3 = Inf;
else
    % find the propagation direction n=V^-T*(d-d3*1)
    v_nan = isnan(V);
    if (v_nan(1,1) | v_nan(1,2) | v_nan(2,1) | v_nan(2,2))
    else
      n=(V^(-1))'*(d-d_3*[1;1]);
      if (Q*V'*n<0) % monotonically condition violated,do dijkstra
          d_x1 = d(1) + normal(1)*Velocity(1);
          d_x2 = d(2) + normal(2)*Velocity(2);
          d_3 = min(d_x1,d_x2);
      end
      d_3 = d_3/Velocity(3); % normalize for output velocity
    end
end

end

