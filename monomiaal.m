function [ M ] = monomiaal( x,n )
M=zeros(length(x),n+1);
for i=1:length(x)
    for j=1:n+1
        M(i,j)=x(i)^(j-1);
    end
end
end

