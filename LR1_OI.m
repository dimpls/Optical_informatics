clear all

a = 1;
b = 10;
c = 5;
p = 0;
q = 3;
n = 1000;
m = 1000;
alfa = 1/1;

hx = (b - a)/n;
x = a : hx : b-hx/2;

function fx_result = f_x(x)
    bet = 1/10;
    fx_result = exp(1i * bet * x);
endfunction

func_res = f_x(x);

function draw_graph(x, f, title_name)
  figure();
  plot(x, f);
  xlabel("x");
  ylabel(title_name);
endfunction

ksi_h = (q-p)/m;
for l=1:1:m
  ksi_l(l) = p + l * ksi_h;
end

function core_result = core(alfa, ksi, x)
    core_result = 1i * x^(alfa * ksi - 1);
endfunction

A = zeros(m,n);
for l=1:1:m
  for k=1:1:n
    A(l,k) = core(alfa, ksi_l(l), x(k));
  endfor
endfor


res_integral = A * func_res' * hx;

#draw_graph(x, abs(func_res), "Амплитуда входного сигнала");
#draw_graph(x, arg(func_res), "Фаза входного сигнала");
draw_graph(ksi_l, abs(res_integral), "Амплитуда выходного сигнала");
#draw_graph(ksi_l, arg(res_integral), "Фаза выходного сигнала");


