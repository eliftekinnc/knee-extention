% Extended Kalman Filter for Knee Extension (With Gravity, Torque Estimation)
% State: [theta; omega; torque] (angle, angular velocity, torque)
% Measurements: [a_t; a_n; omega_gyro]
clear all; close all; clc;
load('knee_data.mat')

% Parametreler
dt = 0.0065; % Örnekleme aralığı (s), ~154 Hz
l = 0.4;     % IMU'nun diz ekleminden uzaklığı (m), ölç!
g = 9.81;    % Yerçekimi (m/s^2)
n_steps = 4967; % 33 saniye, ~154 Hz için (verine göre ayarla)
n = 4967;

% Veri hazırlığı
w = (-1) * knee.GYROXdegs(1:n); % deg/s, negatif işaretli
w_filt = medfilt1(w, 15);

a_n = (-1) * knee.ACCYG(1:n); % Normal ivme, negatif işaretli
a_n_filt = medfilt1(a_n, 15);

a_t = knee.ACCZG(1:n); % Tanjansiyel ivme, emin değilsen kontrol et
a_t_filt = medfilt1(a_t, 15);

% ilk 1 saniyede çıkar
w_filt = w_filt - mean(w_filt(1:150));
a_n_filt = a_n_filt - mean(a_n_filt(1:150));
a_t_filt = a_t_filt - mean(a_t_filt(1:150));

% Custom Linear Matrices (EKF için F lineer, ama H nonlinear)
F = [1, 0.0065, 0;    % theta = theta + omega * dt
     0, 1, dt;        % omega = omega + torque * dt (basit model)
     0, 0, 1];        % torque sabit kalır (noise ile)
B = [0; 1; 0];        % Kontrol matrisi, torque etkisi (u = torque)

% Process ve measurement noise kovaryans
Q = diag([0.0001, 0.01, 0.001]); % [theta, omega, torque]
R = diag([1, 1, 0.01]);          % [a_t, a_n, omega_gyro]

% Başlangıç durumu ve kovaryans
x = [0; 0; 0]; % [theta; omega; torque]
P = eye(3) * 1; % Başlangıç belirsizliği

% Sonuçları saklamak için
estimates = zeros(n_steps, 3); % [theta, omega, torque] tahminleri
innovations = zeros(n_steps, 3); % Innovation için

% EKF Döngüsü
for k = 1:n_steps
    % Predict Adımı (Lineer F ile)
    u = 0; % Tork ölçümü yok, sıfır kabul ediliyor
    x_pred = F * x + B * u; % Lineer tahmin
    P_pred = F * P * F' + Q; % Kovaryans tahmini
    
    % Ölçüm modeli (yerçekimi dahil, nonlinear)
    theta_pred = x_pred(1);
    omega_pred = x_pred(2);
    y_meas = [a_t_filt(k); a_n_filt(k); w_filt(k)]; % Filtreli veriler
    y_pred = [l * x_pred(3) - g * sin(theta_pred); ... % a_t = l * torque - g sin(θ)
              -l * omega_pred^2 - g * cos(theta_pred); ... % a_n = -l ω² - g cos(θ)
              omega_pred]; % omega_gyro = omega
    
    % Jacobian H (nonlinear ölçüm için türev)
    H = [-g * cos(theta_pred), 0, l; ... % ∂a_t/∂θ, ∂a_t/∂ω, ∂a_t/∂τ
         g * sin(theta_pred), -2 * l * omega_pred, 0; ... % ∂a_n/∂θ, ∂a_n/∂ω, ∂a_n/∂τ
         0, 1, 0]; % ∂omega_gyro/∂θ, ∂omega_gyro/∂ω, ∂omega_gyro/∂τ
    
    % Update Adımı (EKF)
    y_hat = y_meas - y_pred; % Innovation
    S = H * P_pred * H' + R; % Innovation kovaryans
    K = P_pred * H' / (S + eps); % Kalman kazancı, numerik stabilite için eps
    x = x_pred + K * y_hat; % Update
    P = (eye(3) - K * H) * P_pred;
    
    % Sonuçları sakla
    estimates(k, :) = x';
    innovations(k, :) = y_hat';
end

% Sonuçları çiz
t = (0:n_steps-1) * dt;
figure;
subplot(3, 1, 1);
plot(t, estimates(:, 1), 'b', 'LineWidth', 2);
title('Tahmin Edilen Açı (theta, rad)');
xlabel('Zaman (s)'); ylabel('Açı (rad)');

subplot(3, 1, 2);
plot(t, w, 'r--', 'LineWidth', 1, 'DisplayName', 'Jiroskop Verisi');
hold on;
plot(t, estimates(:, 2), 'b', 'LineWidth', 2, 'DisplayName', 'Tahmin (EKF)');
title('Açısal Hız (omega, rad/s)');
xlabel('Zaman (s)'); ylabel('Açısal Hız (rad/s)');
legend;

subplot(3, 1, 3);
plot(t, estimates(:, 3), 'b', 'LineWidth', 2);
title('Tahmin Edilen Tork (torque)');
xlabel('Zaman (s)'); ylabel('Tork');

% Innovation grafiği
figure;
subplot(3, 1, 1);
plot(t, innovations(:, 1), 'k', 'LineWidth', 1);
title('Innovation (a_t hatası)');
xlabel('Zaman (s)'); ylabel('Hata (m/s^2)');

subplot(3, 1, 2);
plot(t, innovations(:, 2), 'k', 'LineWidth', 1);
title('Innovation (a_n hatası)');
xlabel('Zaman (s)'); ylabel('Hata (m/s^2)');

subplot(3, 1, 3);
plot(t, innovations(:, 3), 'k', 'LineWidth', 1);
title('Innovation (omega_gyro hatası)');
xlabel('Zaman (s)'); ylabel('Hata (rad/s)');

% Çıktıyı kaydet
save('ekf_with_gravity_estimates.mat', 'estimates', 'innovations');