function [ y ] = evalueer_lagrange( x, f, z )

% initialize y
y = 0;

% iterate through all points
for i = 1:length(x)
    t = f(i);
    
    for j = 1:length(x)
        
        if i ~= j
           t = t .* (z - x(j))/(x(i)-x(j));
        end
    end
    
    y = y + t;
end    

end

