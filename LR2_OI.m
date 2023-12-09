close all;
clear all;

% Задание 1
function F = bpf(f, N, M, h_x)
  zeros_to_add = zeros(1, (M - N) / 2);
  f = horzcat(zeros_to_add, f, zeros_to_add);

  f = fftshift(f);
  F = fft(f) * h_x;
  F = fftshift(F);

  F = F(length(F) / 2 - N / 2 + 1 : length(F) / 2 + N / 2);
end
#


function F = numericFT(f, x_area, u_area, hx, core)
    F = zeros(1, length(u_area));

    for k = 1:length(u_area)
        sum = 0;
        for n = 1:length(x_area)
            sum += f(n) * core(x_area(n), u_area(k));
        end
        F(k) = sum * hx;
    end
end
%
function F = bpf_dd(f, N, M, h_x)
  F = zeros(N, N);
  for i = 1 : N
    F(i, :) = bpf(f(i, :), N, M, h_x);
  endfor
  F = F.';
  for j = 1 : N
    F(j, :) = bpf(F(j, :), N, M, h_x);
  endfor
  F = F.'
end

function y = analitic_fun(x)
    y = zeros(size(x));
    y = 2 ./ (1 + 4*pi^2*x.^2)
end


N = pow2(6);
M = pow2(8);

a = 5;
h_x = 2 * a / N;

x_area = -a : h_x : a - h_x;
b_wave = (N .^ 2) / (4 * a * M);
h_u = (2 * b_wave) / N;
u_area = -b_wave : h_u : b_wave - h_u;

fun_lam = exp(-(x_area .^ 2));

F_gauss = bpf(fun_lam, N, M, h_x);


% Задание 2-3
% Амплитуда гауссова пучка
figure;
plot(x_area, abs(fun_lam), 'LineWidth', 2, 'Color', 'g');
title('Амплитуда входного сигнала');
xlabel('x');
ylabel('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

% Фаза гауссова пучка
figure;
plot(x_area, angle(fun_lam), 'LineWidth', 2, 'Color', 'm');
title('Фаза входного сигнала');
xlabel('x');
ylabel('Фаза');
grid on;
set(gca, 'FontSize', 12);

% Амплитуда гауссова пучка после БПФ
figure;
plot(u_area, abs(F_gauss), 'LineWidth', 2, 'Color', 'g');
title('Амплитуда выходного сигнала после БПФ');
xlabel('u');
ylabel('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

% Фаза гауссова пучка после БПФ
figure;
plot(u_area, angle(F_gauss), 'LineWidth', 2, 'Color', 'm');
title('Фаза выходного сигнала после БПФ');
xlabel('u');
ylabel('Фаза');
grid on;
set(gca, 'FontSize', 12);

% Задание 4
core = @(x, u) exp(-2 * pi * i * u * x);

F_gauss_a = numericFT(fun_lam, x_area, u_area, h_x, core);

% Амплитуда гауссова пучка после A
figure;
plot(u_area, abs(F_gauss_a), 'LineWidth', 2, 'Color', 'g');
title('Амплитуда гауссова пучка (финитное преобразование Фурье)');
xlabel('u');
ylabel('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

% Фаза гауссова пучка после A
figure;
plot(u_area, angle(F_gauss_a), 'LineWidth', 2, 'Color', 'm');
title('Фаза гауссова пучка (финитное преобразование Фурье)');
xlabel('u');
ylabel('Фаза');
grid on;
set(gca, 'FontSize', 12);

%Задание 5
%Амплитуда

% График Амплитуды
figure('Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
plot(u_area, abs(F_gauss), 'LineWidth', 2, 'Color', 'g');
title('Амплитуда - БПФ');
xlabel('u');
ylabel('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

subplot(1, 2, 2);
plot(u_area, abs(F_gauss_a), 'LineWidth', 2, 'Color', 'r');
title('Амплитуда - Финитное преобразование Фурье');
xlabel('u');
ylabel('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

% График Фазы
figure('Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
plot(u_area, angle(F_gauss), 'LineWidth', 2, 'Color', 'b');
title('Фаза - БПФ');
xlabel('u');
ylabel('Фаза');
grid on;
set(gca, 'FontSize', 12);

subplot(1, 2, 2);
plot(u_area, angle(F_gauss_a), 'LineWidth', 2, 'Color', 'm');
title('Фаза - Финитное преобразование Фурье');
xlabel('u');
ylabel('Фаза');
grid on;
set(gca, 'FontSize', 12);


%Задание 6
f_own = exp(-abs(x_area));

F_analitic = ((1 + 4 * pi * pi * u_area .* u_area) / 2) .^ (-1);

F_gauss_own = bpf(f_own, N, M, h_x);


%Задание 6 - 7
%Амплитуда

figure('Position', [100, 100, 1200, 600]);

subplot(2, 2, 1);
plot(x_area, abs(f_own), 'LineWidth', 2, 'Color', 'm');
title('Амплитуда - своего входного сигнала');
xlabel('x');
ylabel('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

subplot(2, 2, 2);
plot(x_area, angle(f_own), 'LineWidth', 2, 'Color', 'm');
title('Фаза - своего входного сигнала');
xlabel('x');
ylabel('Фаза');
grid on;
set(gca, 'FontSize', 12);

subplot(2, 2, 3);
plot(u_area, abs(F_gauss_own), 'g-', 'LineWidth', 2);
title('Амплитуда - своего входного сигнала после БПФ');
xlabel('u');
ylabel('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

subplot(2, 2, 4);
plot(u_area, angle(F_gauss_own), 'Color', 'm', 'LineWidth', 2);
title('Фаза - своего входного сигнала после БПФ');
xlabel('u');
ylabel('Фаза');
grid on;
set(gca, 'FontSize', 12);



figure('Position', [100, 100, 1200, 600]);

plot(u_area, abs(F_gauss_own), 'g-', 'LineWidth', 2);
hold on;
plot(u_area, abs(F_analitic), 'b--', 'LineWidth', 2);
hold off;

title('Сравнение амплитуд');
xlabel('u');
ylabel('Амплитуда');
grid on;
legend('БПФ', 'Аналитическое', 'Location', 'best');
set(gca, 'FontSize', 12);

%Фаза
figure('Position', [100, 100, 1200, 600]);

plot(u_area, angle(F_gauss_own), 'g-', 'LineWidth', 2);
hold on;
plot(u_area, angle(F_analitic), 'b--', 'LineWidth', 2);
hold off;

title('Сравнение фаз');
xlabel('u');
ylabel('Фаза');
grid on;
legend('БПФ', 'NUM', 'Location', 'best');
set(gca, 'FontSize', 12);


% Задание 8
% Амплитуда
%figure('Position', [100, 100, 1200, 600]);

% Амплитуда
%subplot(2, 1, 1);
%plot(u_area, abs(F_gauss_own), 'LineWidth', 2, 'Color', 'g');
%title('Амплитуда гауссова пучка аналитическое');
%xlabel('u');
%ylabel('Амплитуда');
%grid on;
%set(gca, 'FontSize', 12);

% Фаза
%subplot(2, 1, 2);
%plot(u_area, angle(F_gauss_own), 'LineWidth', 2, 'Color', 'm');
%title('Фаза гауссова пучка аналитическое');
%xlabel('u');
%ylabel('Фаза');
%grid on;
%set(gca, 'FontSize', 12);

% 2D

%Задание 9

N = pow2(6);
M = pow2(8);

a = 5;

h_x = 2 * a / N;

x_area = -a : h_x : a - h_x;
f_lamd_dd = @(x, y) exp(-abs(x)) * exp(-abs(y));
fun_gauss_lam = @(x, y) exp(-(x.^2)-(y.^2));

f_dd = zeros(N, N);
f_gauss = zeros(N, N);

for i =  1 : N
  for j = 1 : N
    f_dd(i,j) = f_lamd_dd(x_area(i), x_area(j));
    f_gauss(i,j) = fun_gauss_lam(x_area(i), x_area(j));
  endfor
endfor


F_DD = bpf_dd(f_dd, N, M, h_x);
F_gauss_DD = bpf_dd(f_gauss, N, M, h_x);


figure('Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
imagesc(x_area, x_area, abs(f_gauss));
title('Амплитуда входного сигнала (ГП)');
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 12);

subplot(1, 2, 2);
imagesc(x_area, x_area, angle(f_gauss));
title('Фаза входного сигнала (ГП)');
xlabel('x');
ylabel('y');
colormap hot;
set(gca, 'FontSize', 12);



figure('Position', [100, 100, 1200, 600]);
subplot(1, 2, 1);
imagesc(u_area, u_area, abs(F_gauss_DD));
title('Амплитуда ГП после БПФ');
xlabel('x');
ylabel('y');
colormap hot;
set(gca, 'FontSize', 12);

subplot(1, 2, 2);
imagesc(u_area, u_area, angle(F_gauss_DD));
title('Фаза ГП после БПФ');
xlabel('x');
ylabel('y');
colormap hot;
set(gca, 'FontSize', 12);


figure('Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
imagesc(x_area, x_area, abs(f_dd));
title('Амплитуда входного сигнала exp(-|x|)exp(-|y|)');
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 12);

% Второй график амплитуды
subplot(1, 2, 2);
imagesc(x_area, x_area, angle(f_dd));
title('Фаза входного сигнала exp(-|x|)exp(-|y|)');
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 12);
colormap hot;



figure('Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
imagesc(u_area, u_area, abs(F_DD));
title('Амплитуда входного сигнала exp(-|x|)exp(-|y|) после БПФ');
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 12);

% Второй график амплитуды
subplot(1, 2, 2);
imagesc(u_area, u_area, angle(F_DD));
title('Фаза входного сигнала exp(-|x|)exp(-|y|) после БПФ');
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 12);
colormap hot;


analytical_2D = zeros(N, N);

for i = 1:N
    for j = 1:N
        analytical_2D(i, j) = F_analitic(i) * F_analitic(j);
    end
end
analytical_2D
figure('Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
imagesc(u_area, u_area, abs(analytical_2D));
title('Амплитуда exp(-|x|)exp(-|y|) аналитически');
xlabel('x');
ylabel('y');
colormap hot;
set(gca, 'FontSize', 12);

subplot(1, 2, 2);
imagesc(u_area, u_area, angle(analytical_2D));
title('Фаза exp(-|x|)exp(-|y|) аналитически');
xlabel('x');
ylabel('y');
colormap hot;
set(gca, 'FontSize', 12);




















