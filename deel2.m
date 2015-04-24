%% 3.1 Equidistante punten en het Runge fenomeen

clc
f = @fx;

nbpoints = 10000;
z = linspace(-1,1,nbpoints);

y = fx(z);
% TODO
x = linspace(-1,1,9);
polynom8 = evalueer_lagrange(x, f, z);
x = linspace(-1,1,11);
polynom10 = evalueer_lagrange(x, f, z);
x = linspace(-1,1,13);
polynom12 = evalueer_lagrange(x, f, z);
x = linspace(-1,1,15);
polynom14 = evalueer_lagrange(x, f, z);

figure(1),clf
plot(z,y)
hold on
%plot(z, polynom8)
plot(z, polynom10)
%plot(z, polynom12)
%plot(z, polynom14)
legend('f(x)','L8', 'L10', 'L12', 'L14')
xlabel('x-waarden in [-1,1]') % x-axis label
ylabel('Veelterminterpolaties') % y-axis label
