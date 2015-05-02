% methode van Newton-Raphson: exact
x_ster=x_diff(2,end);
x_start=0.2;
n=0.05;
t=750;
p_n=[1:t];
for k=1:t
    p_n(k)=abs((x_diff(5,k+1)-x_ster)/((x_diff(5,k)-x_ster)^n));
end

p_n