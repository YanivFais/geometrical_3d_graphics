function [ x1,x2] = quadratic( a,b,c )
%compute [x1,x2] = a*x^2 + b*x+c where all variables are scalars

if (a==0)
    x1 = -c/b;
    x2 = x1;
else
    x1 = -(b + (b^2 - 4*a*c)^(1/2))/(2*a);
    x2 = -(b - (b^2 - 4*a*c)^(1/2))/(2*a);
end

end

