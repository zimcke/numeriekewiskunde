function [ y ] = fx( x )
y = zeros(length(x), 1);
for i=1:length(x)
    y(i) = 1 /(1 + (25 * (x(i)^2)));
end
end

