clc; clear; close all;
mode = input('Select Signal Type (1 = Discrete, 2 = Continuous): ');

if mode == 1
    x = input('Enter discrete x[n] (e.g., [1 2 3]): ');
    h = input('Enter discrete h[n] (e.g., [1 1]): ');
    dt = 1; tx = 0:length(x)-1; th = 0:length(h)-1;
else
    t = input('Enter time vector t (e.g., -1:0.05:4): ');
    dt = t(2) - t(1);
    fx = str2func(['@(t) ' input('Enter x(t) (e.g., double(t>=0 & t<=2)): ', 's')]);
    fh = str2func(['@(t) ' input('Enter h(t) (e.g., exp(-t).*(t>=0)): ', 's')]);
    x = fx(t); h = fh(t); tx = t; th = t;
end

Nx = length(x); Nh = length(h); Ny = Nx + Nh - 1;
x_pad = [zeros(1, Nh-1), x, zeros(1, Nh-1)];
t_y   = (tx(1) + th(1)) + (0:Ny-1)*dt;
tau   = (tx(1) - (Nh-1)*dt) + (0:length(x_pad)-1)*dt;
y = zeros(1, Ny); y_full = conv(x, h)*dt; h_flip = fliplr(h);

pb = max(abs(x))*max(abs(h)); if pb == 0, pb = 1; end
p_ylim = [-0.15 - (min(x)<0 || min(h)<0), 1.15] * pb;

figure('Name', 'Universal Convolution Animator', 'Color', 'w');
step = 1 + 2*(mode == 2); prev_i = 0;

for i = unique([1:step:Ny, Ny])
    h_shifted = [zeros(1, i-1), h_flip, zeros(1, Ny-i)];
    prod_sig  = x_pad .* h_shifted;
    for j = prev_i+1:i
        y(j) = sum(x_pad .* [zeros(1, j-1), h_flip, zeros(1, Ny-j)]) * dt;
    end
    prev_i = i;

    subplot(3,1,1);
    if mode == 1, stem(tau, x_pad, 'b', 'filled'); hold on; stem(tau, h_shifted, 'r--'); hold off;
    else, plot(tau, x_pad, 'b', tau, h_shifted, 'r--', 'LineWidth', 1.5); end
    title(sprintf('Flip & Shift at t = %.2f', t_y(i)));
    legend('x(\tau)', 'h(t - \tau)'); grid on; xlim([tau(1) tau(end)]);

    subplot(3,1,2);
    if mode == 1
        stem(tau, prod_sig, 'm', 'filled', 'LineWidth', 1.5);
        title(sprintf('Multiply: x[k] \\cdot h[n - k]   |   Current Sum = %.2f', y(i)));
    else
        area(tau, prod_sig, 'FaceColor', [.8 .2 .8], 'FaceAlpha', .35, 'EdgeColor', [.5 0 .5], 'LineWidth', 1.5);
        title(sprintf('Multiply: x(\\tau) \\cdot h(t - \\tau)   |   Shaded Area = %.2f', y(i)));
    end
    yline(0, 'k-', 'LineWidth', 0.8); grid on; xlim([tau(1) tau(end)]); ylim(p_ylim);
    xlabel('\tau'); ylabel('Product');

    subplot(3,1,3);
    if mode == 1, stem(t_y(1:i), y(1:i), 'g', 'filled', 'LineWidth', 1.5);
    else, plot(t_y(1:i), y(1:i), 'g', 'LineWidth', 2); end
    title(sprintf('Sum / Integrate: Growing Output y = %.2f', y(i)));
    grid on; xlim([t_y(1) t_y(end)]); ylim([min(y_full)-0.5, max(y_full)+0.5]);

    drawnow; pause(0.25);
end
