function [ T, dT ] = chebychev2( x,n )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


if(n == 0)
            T(i,j) = 1;
            dT(i,j) = 0;
end
if(n == 1)
            T(i,j) = x(i);
            dT(i,j) = 1;


end

end

