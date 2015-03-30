function [ T, dT ] = chebyshev( x,n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
T = double(zeros(length(x),n))
dT = double(zeros(length(x),n))

if(n == 0)
    for j = drange(0:n)
        for i = drange(1:size(x))
            T(i,j) = 1;
            dT(i,j) = 0;
        end
    end
end
if(n == 1)
     for j = drange(n:0)
        for i = drange(size(x):1)
            T(i,j) = x(i);
            dT(i,j) = 1;
        end
     end
else
    for j = drange(n:0)
        for i = drange(size(x):1)
               T(i,j) = 2*x*chebyshev(x,j-1) - chebyshev(x,j-2);
               dT(i,j) = 1;
        end
    end

end
end

