%yerçekimsiz ortamda çembersel hareket state'e tork eklendi + ivmeölçer
%verileir eklendi


clear all; close all; clc;
load('knee_data.mat')

g = 9.81;
dt = 0.0067; 
l = 0.4;
m = 5;
I = 1/3 * m * l^2;
n = 4967;

% Veri hazırlığı
w = (-1) * knee.GYROXdegs(1:n); % deg/s, negatif işaretli
w_filt = medfilt1(w, 15);

a_n = (-1) * g * knee.ACCYG(1:n); % Normal ivme, negatif işaretli (kullanılmıyor)
a_n_filt = medfilt1(a_n, 15);

a_t = g * knee.ACCZG(1:n); % Tanjansiyel ivme (m/s^2)
a_t_filt = medfilt1(a_t, 15);

% İlk 1 saniyede ortalama çıkar (bias kaldırma)
w_filt = w_filt - mean(w_filt(1:150));
a_n_filt = a_n_filt - mean(a_n_filt(1:150));
a_t_filt = a_t_filt - mean(a_t_filt(1:150));

% Sistem Matrisleri (Lineer)
F = [1, dt, 0; 
     0, 1, dt/I; 
     0, 0, 1];
B = [0; 0; 0]; % Kontrol inputu yok

% Ölçüm Matrisi (omega ve alpha = tork / I)
H = [0, 1, 0;    % omega ölçümü
     0, 0, l/I]; % alpha = tork / I (a_t / l = alpha)

% Process ve Measurement Noise Kovaryansları
Q = diag([0.01, 0.01, 0.01]); % [theta, omega, tork] için process noise
R = diag([0.1, 0.1]);         % [omega_gyro (rad/s), alpha_accel (rad/s^2)] için measurement noise
% R değerlerini verilerin gürültüsüne göre ayarlayabilirsiniz. Burada örnek olarak 0.1 kullandım.

% Başlangıç Durumu ve Kovaryans
x = [0; 0; 0]; % [theta (rad); omega (rad/s); tork (Nm)]
P = eye(3) * 1; % Başlangıç belirsizliği

% Sonuçları Saklamak İçin
estimates = zeros(n, 3); % [theta, omega, tork] tahminleri
innovations = zeros(n, 2); % [omega_hat, alpha_hat] innovation'lar

% Linear Kalman Filter Döngüsü
for k = 1:n
    % Predict Adımı
    x_pred = F * x;
    P_pred = F * P * F' + Q;
    
    % Update Adımı
    a_t_meas = a_t_filt(k) ;
    y_meas = [deg2rad(w_filt(k)); a_t_meas];
    y_pred = H * x_pred;
    y_hat = y_meas - y_pred;
    S = H * P_pred * H' + R;
    K = P_pred * H' / (S + eps * eye(size(S))); % Kalman kazancı, stabilite için eps
    x = x_pred + K * y_hat;
    P = (eye(3) - K * H) * P_pred;
    
    % Sonuçları Sakla
    estimates(k, :) = x';
    innovations(k, :) = y_hat';
end

t = (0:n-1) * dt;

% Plot Tahminler
figure;
subplot(3, 1, 1);
plot(t, estimates(:, 1), 'b', 'LineWidth', 2);
grid on;
title('Tahmin Edilen Açı \theta (rad)');
xlabel('Zaman (s)'); ylabel('\theta (rad)');

subplot(3, 1, 2);
plot(t, deg2rad(w_filt), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Jiroskop (rad/s)'); hold on;
plot(t, estimates(:, 2), 'b', 'LineWidth', 2, 'DisplayName', 'KF \omega (rad/s)');
grid on;
title('Açısal Hız \omega (rad/s)');
xlabel('Zaman (s)'); ylabel('\omega (rad/s)');
legend('Location', 'best');

subplot(3, 1, 3);
plot(t, estimates(:, 3), 'b', 'LineWidth', 2);
grid on;
title('Tahmin Edilen Tork (Nm)');
xlabel('Zaman (s)'); ylabel('Tork (Nm)');

% Plot Innovation'lar
figure;
subplot(2, 1, 1);
plot(t, innovations(:, 1), 'b', 'LineWidth', 1);
grid on;
title('Innovation (Gyro \omega)');
xlabel('Zaman (s)'); ylabel('Hata (rad/s)');

subplot(2, 1, 2);
plot(t, innovations(:, 2), 'b', 'LineWidth', 1);
grid on;
title('Innovation (Accel \alpha = a_t / l)');
xlabel('Zaman (s)'); ylabel('Hata (rad/s^2)');

save('kf_linear_gyro_accel.mat', 'estimates', 'innovations', 't');