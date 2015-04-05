function [ T , dT ] = chebyshev( x,n )
T=zeros(length(x),n+1);
dT=zeros(length(x),n+1);
for i=1:length(x)
    for j =1:n+1
        if j==1
            T(i,j)=1;
            dT(i,j)=0;
        elseif j==2;
            T(i,j)=x(i);
            dT(i,j)=1;
        else
            T(i,j)=2*x(i)*T(i,j-1)-T(i,j-2);
            dT(i,j)=2*x(i)*dT(i,j-1)+2*T(i,j-1)-dT(i,j-2);
        end
    end
end
end

