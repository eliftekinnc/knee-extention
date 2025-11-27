clear all
load('knee_data.mat')

w = (-1)*knee.GYROXdegs(1:4967); %deg/s
z = deg2rad(w);  % deg/s -> rad/s
z_static = z(1:15);   % örnek: ilk ~5 sn, dt=0.0065 ise
R = var(z_static);                  % (rad/s)^2


dt = 0.0065;                 % örnekleme periyodu (s)
alpha_std = 10;              % tipik açısal ivme (rad/s^2) ~ kaba tahmin
Q = (alpha_std*dt)^2;        % süreç gürültüsü önerisi
R = var(z_static);           % statik durumda toplanan gyro varyansı

[omega_hat, K, P] = kalman1d_gyro(z, Q, R);

t = (0:numel(z)-1)'*dt;   % zaman ekseni
plot(t, z, 'r--', t, omega_hat, 'b', 'LineWidth',1.2)
legend('Ham gyro (z)', 'Kalman çıktı (\omegâ)')
xlabel('zaman (s)'); ylabel('rad/s'); grid on
