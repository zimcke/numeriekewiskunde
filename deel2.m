%% 3.1 Equidistante punten en het Runge fenomeen

clc
f = @fx;

nbpoints = 10000;
x = linspace(-1,1,nbpoints);

y = fx(x);
% TODO
%polynom8 = 
%polynom10 = 
%polynom12 = 
%polynom14 = 

figure(1),clf
plot(x,y)
hold on
plot(x, polynom8)
plot(x, polynom10)
plot(x, polynom12)
plot(x, polynom14)
legend('f(x)','L8', 'L10', 'L12', 'L14')
xlabel('x-waarden in [-1,1]') % x-axis label
ylabel('Veelterminterpolaties') % y-axis label
