% input x = the horizontal points in which the monomial polynomial are to be evaluated
% input n = the maximum degree to which to compute the interpolants
% output M = the monomial interpolants in all x values, to degree n
function [ M ] = monomiaal( x,n )

% initialize M
M = zeros(length(x),n+1);

% iterate through horizontal points
for i = 1:length(x)
    
    % iterate through all degrees
    for j = 1:n+1
        M(i,j)=x(i)^(j-1);
    end
end
end

