function [ nodes ] = chebyshev_nodes( i )
% returns a vector of length i containing Chebyshev nodes
nodes = cos( (2*(1:i)-1)*pi/2/(i));
end

