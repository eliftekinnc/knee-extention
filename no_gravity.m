%yerçekimsiz ortamda çembersel hareket
clear all; close all; clc;
load('knee_data.mat')

g = 9.81;
dt = 0.0067; 
l = 0.4;
n = 4967;

% Veri hazırlığı
w = (-1) * knee.GYROXdegs(1:n); % deg/s, negatif işaretli
w_filt = medfilt1(w, 15);

a_n = (-1) * g * knee.ACCYG(1:n); % Normal ivme, negatif işaretli
a_n_filt = medfilt1(a_n, 15);

a_t = g * knee.ACCZG(1:n); % Tanjansiyel ivme, emin değilsen kontrol et
a_t_filt = medfilt1(a_t, 15);

% ilk 1 saniyede çıkar
w_filt = w_filt - mean(w_filt(1:150));
a_n_filt = a_n_filt - mean(a_n_filt(1:150));
a_t_filt = a_t_filt - mean(a_t_filt(1:150));


% Custom Linear Matrices
F = [1, dt; 0, 1];
B = [0; dt/l]; % kontrol inputu a_n
H = [0, 1]; % ölçüm omega

% Process ve measurement noise kovaryans (lineer için)
Q = diag([0.01, 0.01]); % [theta, omega]
R = 0.1  % [omega_gyro]

% Başlangıç durumu ve kovaryans
x = [0; 0]; % [theta; omega]
P = eye(2) * 1; % Başlangıç belirsizliği

% Sonuçları saklamak için
estimates = zeros(n, 2); % [theta, omega] tahminleri
innovations = zeros(n,1);

% Linear Kalman Filter Döngüsü
for k = 1:n
    % Predict Adımı 
    x_pred = F * x; % Lineer tahmin
    P_pred = F*P*F' + Q ;  % tahmin kovaryansı
    
    
    % Update Adımı (Lineer)
    y_meas = deg2rad(w_filt(k));
    y_pred = H * x_pred;
    y_hat = y_meas - y_pred;
    S = H * P_pred * H' + R; % Sadece gyro ölçümü için (H [0, 1])
    K = P_pred * H' / (S + eps); % Kalman kazancı, numerik stabilite için eps
    x = x_pred + K * y_hat; % durum güncellemesi
    P = (eye(2) - K * H) * P_pred; 
    
    % Sonuçları sakla
    estimates(k, :) = x';
    innovations(k, :) = y_hat';
end

t = (0:n-1) * dt;

figure;
subplot(2, 1, 1);
plot(t, estimates(:, 1), 'b', 'LineWidth', 2);
grid on;
title('Tahmin Edilen Açı \theta (rad)');
xlabel('Zaman (s)'); ylabel('\theta (rad)');

subplot(2, 1, 2);
plot(t, deg2rad(w_filt), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Jiroskop (rad/s)'); hold on;
plot(t, estimates(:, 2), 'b', 'LineWidth', 2, 'DisplayName', 'KF \omega (rad/s)');
grid on;
title('Açısal Hız \omega (rad/s)');
xlabel('Zaman (s)'); ylabel('\omega (rad/s)');
legend('Location','best');




figure;
plot(t, innovations, 'k', 'LineWidth', 1);
grid on;
title('Innovation (gyro \omega)');
xlabel('Zaman (s)'); ylabel('Hata (rad/s)');


save('kf_linear_gyro_only.mat', 'estimates', 'innovations', 't');
