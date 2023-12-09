close all;
clear all;

R = 5;
m = 2;
n = 256;
N = 2 * n + 1;
h_r = R / n;
r_h = 0 : h_r : (R - h_r / 2);
M = 4096;
b = (n ^ 2) / (4 * R * M);
h_x = (2 * R) / n;

function res = func_ro(ro)
    res = (5 * (ro .^ 5) - 4 * ro .^ 3);
end

values_func_ro = func_ro(r_h);

%Задание 1

figure('Position', [100, 100, 1200, 600]);
plot(r_h, abs(values_func_ro), 'LineWidth', 2, 'Color', 'g');
title('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

figure('Position', [100, 100, 1200, 600]);
plot(r_h, angle(values_func_ro), 'LineWidth', 2, 'Color', 'g');
title('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

%Задание 2

function a = two_d_array_of_function_values(values, N, m)
    n = N - 1;
    a = zeros(2 * N + 1, 2 * N + 1);
    for j = 1:2 * N + 1
        for k = 1:2 * N + 1
            alpha = round(sqrt((j - n - 1) ^ 2 + (k - n - 1) ^ 2));
            if (alpha <= n)
                a(j, k) = values(alpha + 1) * exp(1i * m * atan2(k - n - 1, j - n - 1));
            end
        end
    end
end

two_d_array_of_func = two_d_array_of_function_values(values_func_ro, n, m);

function plot_2d_function(input_function, R, amplitude_title, phase_title)
    fig = figure;
    ax1 = subplot(1, 2, 1);
    ax2 = subplot(1, 2, 2);

    imagesc(ax1, -R:R, -R:R, abs(input_function));
    colormap(ax1, 'Winter'); %Autumn
    title(ax1, amplitude_title);
    colorbar(ax1);

    imagesc(ax2, -R:R, -R:R, angle(input_function));
    colormap(ax2, 'Winter');
    title(ax2, phase_title);
    colorbar(ax2);
end

%plot_2d_function(two_d_array_of_func, R, 'Амплитуда', 'Фаза');
figure('Position', [100, 100, 600, 600]);
imagesc(-R:R, -R:R, abs(two_d_array_of_func));
title('Амплитуда');
colormap('Winter');
colorbar;

figure('Position', [100, 100, 600, 600]);
imagesc(-R:R, -R:R, angle(two_d_array_of_func));
title('Фаза');
colormap('Winter');
colorbar;

%Задание 3

function integral = hankel_transformation(input_function, m, input_area, h_r)
    x = input_area;
    n = length(x);
    integral = complex(zeros(1, n));
    for i = 1:n
        integral(i) = sum(input_function .* besselj(m, 2 * pi * x * x(i)) .* x * h_r);
    end
    integral = integral * (2 * pi / 1i^m);
end

tic;
H = hankel_transformation(values_func_ro, m, r_h, h_r);
elapsedTimeH = toc

figure('Position', [100, 100, 600, 600]);
plot(r_h, abs(H), 'LineWidth', 2, 'Color', 'g');
title('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

figure('Position', [100, 100, 600, 600]);
plot(r_h, angle(H), 'LineWidth', 2, 'Color', 'g');
title('Фаза');
grid on;
set(gca, 'FontSize', 12);

H_2D_array = two_d_array_of_function_values(H, n, m);

%plot_2d_function(H_2D_array, R, 'Амплитуда', 'Фаза');

figure('Position', [100, 100, 600, 600]);
imagesc(-R:R, -R:R, abs(H_2D_array));
title('Амплитуда');
colormap('Winter');
colorbar;

figure('Position', [100, 100, 600, 600]);
imagesc(-R:R, -R:R, angle(H_2D_array));
title('Фаза');
colormap('Winter');
colorbar;
%Задание 4

function F = bpf(f, N, M, h_x)
  zeros_to_add = zeros(1, round((M - N) / 2));
  f = horzcat(zeros_to_add, f, zeros_to_add);

  f = fftshift(f);
  F = fft(f) * h_x;
  F = fftshift(F);

  F = F(length(F) / 2 - N / 2 + 1 : length(F) / 2 + N / 2);
end

function F = bpf_dd(f, N, M, h_x)
  F = zeros(N, N);
  for i = 1 : N
    F(i, :) = bpf(f(i, :), N, M, h_x);
  endfor
  F = F.';
  for j = 1 : N
    F(j, :) = bpf(F(j, :), N, M, h_x);
  endfor
  F = F.';
end

tic;
FFT_2D = bpf_dd(two_d_array_of_func, N, M, h_x);
elapsedTimeBPF = toc

%plot_2d_function(FFT_2D, b, 'Амплитуда', 'Фаза');

h_u = (2 * b) / N;
u_area = -b : h_u : b - h_u;

figure('Position', [100, 100, 600, 600]);
imagesc(u_area, u_area, abs(FFT_2D));
title('Амплитуда');
colormap('Winter');
colorbar;

figure('Position', [100, 100, 600, 600]);
imagesc(u_area, u_area, angle(FFT_2D));
title('Фаза');
colormap('Winter');
colorbar;

