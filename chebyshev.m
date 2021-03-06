% input x = the horizontal points in which the Chebychev polynomial are to be evaluated
% input n = the maximum degree to which to compute the interpolants
% output T = the Chebyshev interpolants in all x values, to degree n
% output dT = derivative of the Chebychev interpolants in all x values, to
% degree n
function [ T , dT ] = chebyshev( x , n )
T = zeros(length(x), n+1);
dT = zeros(length(x), n+1);

% iterate through all x values
for i = 1:length(x)
    
    % iteratre through all degrees
    for j = 1:n+1
        
        % if degree == 0
        if j==1
            T(i,j)=1;
            dT(i,j)=0;
            
        % if degree == 2    
        elseif j==2;
            T(i,j)=x(i);
            dT(i,j)=1;
               
        else
            % recursion formula
            T(i,j)=2*x(i)*T(i,j-1)-T(i,j-2);
            dT(i,j)=2*x(i)*dT(i,j-1)+2*T(i,j-1)-dT(i,j-2);
        end
    end
end
end

