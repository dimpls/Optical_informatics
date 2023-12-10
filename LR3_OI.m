close all;
clear all;

R = 5;
m = 2;

n = 256;
N = 2 * n + 1;
M = 4096;

step = R / n;
r_linspace = 0 : step : (R - step / 2);

b_wave = (n ^ 2) / (4 * R * M);
h_x = 2 * R / n;

range = -R:R;

h_u = (2 * b_wave) / N;
u_area = -b_wave : h_u : b_wave - h_u;

function res = func_ro(ro)
    res = (5 * (ro .^ 5) - 4 * ro .^ 3);
end

values_func_ro = func_ro(r_linspace);

%Задание 1

figure('Position', [100, 100, 1200, 600]);
plot(r_linspace, abs(values_func_ro), 'LineWidth', 2, 'Color', 'g');
title('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

figure('Position', [100, 100, 1200, 600]);
plot(r_linspace, angle(values_func_ro), 'LineWidth', 2, 'Color', 'g');
title('Фаза');
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

figure('Position', [100, 100, 600, 600]);
imagesc(range, range, abs(two_d_array_of_func));
title('Амплитуда');
colormap('Winter');
colorbar;

figure('Position', [100, 100, 600, 600]);
imagesc(range, range, angle(two_d_array_of_func));
title('Фаза');
colormap('Winter');
colorbar;

%Задание 3

function integral = H_transformation(input_function, m, area, step)
    x = area;
    n = length(x);
    integral = complex(zeros(1, n));
    for i = 1:n
        integral(i) = sum(input_function .* besselj(m, 2 * pi * x * x(i)) .* x * step);
    end
    integral = integral * (2 * pi / 1i^m);
end

tic;
H = H_transformation(values_func_ro, m, r_linspace, step);
TimeH = toc

figure('Position', [100, 100, 600, 600]);
plot(r_linspace, abs(H), 'LineWidth', 2, 'Color', 'g');
title('Амплитуда');
grid on;
set(gca, 'FontSize', 12);

figure('Position', [100, 100, 600, 600]);
plot(r_linspace, angle(H), 'LineWidth', 2, 'Color', 'g');
title('Фаза');
grid on;
set(gca, 'FontSize', 12);

H_2D_array = two_d_array_of_function_values(H, n, m);


figure('Position', [100, 100, 600, 600]);
imagesc(range, range, abs(H_2D_array));
title('Амплитуда');
colormap('Winter');
colorbar;

figure('Position', [100, 100, 600, 600]);
imagesc(range, range, angle(H_2D_array));
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
TimeBPF = toc

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

