function [ y,dy ] = chebpolyval( c,x )
[T,dT]=chebyshev3(x,length(c)-1);
y=T*c';
y=y';
dy=dT*c';
dy=dy'
end

