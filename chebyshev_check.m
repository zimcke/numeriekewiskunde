% input x = the horizontal points in which the Chebychev polynomial are to be evaluated
% input n = the maximum degree to which to compute the interpolants
% output T = the Chebyshev interpolants in all x values, to degree n
% output dT = derivative of the Chebychev interpolants in all x values, to
% degree n
function [ T , dT ] = chebyshev_check( x , n )
T = zeros(length(x), n+1);
dT = zeros(length(x), n+1);

% iterate through all x values
for i = 1:length(x)
    
    % iteratre through all degrees
    for j = 1:n+1
        
            T(i,j)= cos((j-1)* acos(x(i)));
            dT(i,j)= ((j-1) * sin((j-1) * acos(x(i))))/sin(acos(x(i)));
       
    end
end
end

