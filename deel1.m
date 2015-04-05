
clc
x = linspace(-1,1,10);
c = [1:10]
[T, dT] = chebyshev3(x,5)
[y, dy] = chebpolyval(c,x)