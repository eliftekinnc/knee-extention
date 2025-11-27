% Diz ekstansiyon hareketi için EKF ve animasyon
% Ölçüm: Gerçek sensör verilerinden açısal hız (omega)
% Başlangıç açısı: 0 derece

% Durum geçiş fonksiyonu

function xNext = stateTransitionFcn(x, u, dt)
    % x = [theta; omega], u = kas torku, dt = zaman adımı
    theta = x(1);
    omega = x(2);
    g = 9.81; % Yerçekimi ivmesi (m/s^2)
    L = 0.4; % Bacak uzunluğu (metre)
    m = 5; % Bacak kütlesi (kg)
    I = m * L^2; % Eylemsizlik momenti (kg·m^2)
    xNext = [theta + omega * dt; omega + (u - m*g*L*sin(theta))/I * dt];
end



% Ölçüm fonksiyonu
function y = measurementFcn(x)
    % Sensör, yalnızca açısal hızı (omega) ölçer
    y = x(2);
end


% Jacobian fonksiyonları
function F = stateTransitionJacobianFcn(x, u, dt)
    % Durum geçiş fonksiyonunun Jacobian'ı
    m = 5;
    g = 9.81;
    L = 0.4;
    I = m * L^2;
    theta = x(1);
    F = [1 dt; -m*g*L*cos(theta)/I*dt 1]; % df/dx
end

function H = measurementJacobianFcn(x)
    % Ölçüm fonksiyonunun Jacobian'ı
    % y = omega => dy/dx = [0 1]
    H = [0 1];
end



% Ana script
clear all; close all; clc;
load('knee_data.mat')

% Sensör verilerini yükleme
t = knee.ACCXTimeSeriess(1:4967); % Zaman vektörü
y_meas = (-1) * deg2rad(knee.GYROXdegs(1:4967)); % Derece/s'den radyan/s'ye çevir
N = length(t); % Veri noktası sayısı
dt = t(2) - t(1); % Zaman adımı (sabit olduğu varsayılıyor)

% Parametreler
u = 0.5; % Sabit kas torku (Nm)
Q = diag([0.001 0.001]); % Proses gürültüsü kovaryansı
R = 0.05; % Ölçüm gürültüsü kovaryansı (sensör özelliklerine göre ayarlayın)



% EKF nesnesi oluşturma
initialState = [0; y_meas(1)]; % Başlangıç tahmini: 0 derece, ilk ölçüm
ekf = extendedKalmanFilter(@(x) stateTransitionFcn(x, u, dt), @measurementFcn, initialState);

% EKF özellikleri
ekf.StateCovariance = eye(2) * 0.1; % Başlangıç hata kovaryansı
ekf.ProcessNoise = Q;
ekf.MeasurementNoise = R;
ekf.StateTransitionJacobianFcn = @(x) stateTransitionJacobianFcn(x, u, dt);
ekf.MeasurementJacobianFcn = @measurementJacobianFcn;

% EKF tahminleri
x_est = zeros(2, N); % Tahmin edilen durumlar
x_est(:, 1) = initialState;

% EKF döngüsü
for k = 2:N
    predict(ekf);
    correct(ekf, y_meas(k));
    x_est(:, k) = ekf.State;
end


% Grafiksel sonuçlar
figure('Position', [900, 100, 800, 600]);
subplot(2, 1, 1);
plot(t, x_true(1, :)*180/pi, 'b-', 'DisplayName', 'Gerçek Açı (derece)');
hold on;
plot(t, x_est(1, :)*180/pi, 'r--', 'DisplayName', 'Tahmin Edilen Açı (derece)');
xlabel('Zaman (s)'); ylabel('Açı (derece)');
title('Diz Açısı Tahmini (Başlangıç: 0°, Sensör: Açısal Hız)');
legend;

subplot(2, 1, 2);
plot(t, x_true(2, :)*180/pi, 'b-', 'DisplayName', 'Gerçek Açısal Hız (derece/s)');
hold on;
plot(t, x_est(2, :)*180/pi, 'r--', 'DisplayName', 'Tahmin Edilen Açısal Hız (derece/s)');
plot(t, y_meas*180/pi, 'g:', 'DisplayName', 'Ölçülen Açısal Hız (derece/s)');
xlabel('Zaman (s)'); ylabel('Açısal Hız (derece/s)');
title('Açısal Hız Tahmini ve Sensör Ölçümleri');
legend;

% Hata analizi
angle_error = x_true(1, :) - x_est(1, :);
omega_error = x_true(2, :) - x_est(2, :);
angle_error_cov = sum(angle_error.^2) / N;
omega_error_cov = sum(omega_error.^2) / N;
fprintf('Açı Hata Kovaryansı: %.4f (rad^2)\n', angle_error_cov);
fprintf('Açısal Hız Hata Kovaryansı: %.4f (rad^2/s^2)\n', omega_error_cov);