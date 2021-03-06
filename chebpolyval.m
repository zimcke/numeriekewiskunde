function [ y,dy ] = chebpolyval( c,x )

% corresponding chebychev polynomials
[T,dT]=chebyshev(x,length(c)-1);

% init
y = zeros(length(x), 1)';
dy = zeros(length(x), 1)';

%% p(x)
% iterate through all x values
for i = 1:length(x)
    for j = 1:length(c)
        y(i) = y(i) + c(j).*T(i,j);  
    end
end

%% p'(x)
% iterate through all x values
for i = 1:length(x)
    for j = 1:length(c)
        dy(i) = dy(i) + c(j).*dT(i,j);  
    end
end

% y = T*c';
% y = y';
% dy = dT*c';
% dy = dy'
end

