%% 3.1 Equidistante punten en het Runge fenomeen

clc
f = @fx;

nbpoints = 10000;
z = linspace(-1,1,nbpoints);

y = fx(z);
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
plot(z, polynom8)
plot(z, polynom10)
plot(z, polynom12)
plot(z, polynom14)
legend('f(x)','L8', 'L10', 'L12', 'L14')
xlabel('x-values in [-1,1]') % x-axis label
ylabel('Polynomial interpolation') % y-axis label
print('report\3_1_1','-dpng')

% equidistantial nodes degree 20
x = linspace(-1,1,21);
polynom20 = evalueer_lagrange(x, f, z);
abs_error20 = abs(polynom20 - y);
rel_error20 = abs_error20 ./ y;

% equidistantial nodes degree 50
x = linspace(-1,1,51);
polynom50 = evalueer_lagrange(x, f, z);
abs_error50 = abs(polynom50 - y);
rel_error50 = abs_error50 ./ y;

% chebychev nodes degree 20
xc = 0;
for i = 1:21
   xc(i) = cos(((2*i - 1)*pi)/ (2*21));
end
xc
polynom20c = evalueer_lagrange(xc, f, z);
abs_error20c = abs(polynom20c - y);
rel_error20c = abs_error20c ./ y;

% chebychev nodes degree 50
for i = 1:51
   xc(i) = cos(((2*i - 1)*pi)/ (2*51));
end
polynom50c = evalueer_lagrange(xc, f, z);
abs_error50c = abs(polynom50c - y);
rel_error50c = abs_error50c ./ y;

% plot error
figure(2),clf
semilogy(z, rel_error20)
hold on
semilogy(z, rel_error50)
semilogy(z, rel_error20c)
semilogy(z, rel_error50c)
legend('Equidistant nodes n=20', 'Equidistant nodes n=50', 'Chebychev nodes n=20', 'Chebychev nodes n=50', 'Location','southwest')
xlabel('x-values in [-1,1]') % x-axis label
ylabel('Relative error') % y-axis label
print('report\3_1_2','-dpng')

%% Verschillende basissen

clc

f = @fx;
nbpoints = 10000;
z = linspace(-1,1,nbpoints);
y = fx(z);

xc = 0;
for i = 1:81
   xc(i) = cos(((2*i - 1)*pi)/ (2*51));
end

lagrange80 = evalueer_lagrange(xc, f, z);
monomial80 = polyval(xc, z);
chebpolyval(xc, z);

