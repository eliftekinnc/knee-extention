%% Diz Açısı + Tork EKF (x = [theta; omega; tau])
clear; close all; clc;
load('knee_data.mat');

% --- Fizik parametreleri ---
g     = 9.81;
L     = 0.3;          % (m) etkin moment kolu (kütle merkezi)
L_tam = 0.5;          % (m) segment uzunluğu (varsayım)
m     = 5;            % (kg)
I     = m * L_tam^2;  % basit atalet (istersen 1/3*m*L_tam^2)
b     = 0.0;          % viskoz sönüm (Nms/rad) -> ayarlanabilir
k     = 0.0;          % elastik rijitlik (Nm/rad) -> ayarlanabilir

% --- Veri ---
N   = 4967;
t   = knee.ACCXTimeSeriess(1:N);
w_d = (-1) * knee.GYROXdegs(1:N);      % deg/s (gyro)
a_t = knee.ACCZG(1:N) * g;             % m/s^2 (bilgi amaçlı)
a_n = (-1) * knee.ACCYG(1:N) * g;      % m/s^2 (bilgi amaçlı)

% --- Filtreleme & bias ---
w_filt   = medfilt1(w_d, 15);
a_t_filt = medfilt1(a_t, 15);
a_n_filt = medfilt1(a_n, 15);

w_filt   = w_filt   - mean(w_filt(1:150));
a_t_filt = a_t_filt - mean(a_t_filt(1:150));
a_n_filt = a_n_filt - mean(a_n_filt(1:150));

% --- Radyan ---
w_rad    = deg2rad(w_filt);            % rad/s
theta_int = cumtrapz(t, w_rad);         % rad (referans çizimi)
alpha_kin = (a_t + g*sin(theta_int))/L;     % rad/s^2 (sadece görselleştirme)

% --- Ölçüm ---
y_meas = w_rad;                        % y = omega (rad/s)

% --- Gürültü kovaryansları ---
% (theta, omega, tau) süreç gürültüsü:
q_theta = 1e-4;
q_omega = 1e-3;
q_tau   = 1e-1;  % torkun değişkenliği (↑ daha “oynak” tau)
Q = diag([q_theta, q_omega, q_tau]);

% ölçüm gürültüsü (omega)
R = (0.05)^2;   % (rad/s)^2 yaklaşık

% --- (Opsiyonel) "gerçek" sistem simülasyonu (kıyas için) ---
x_true = zeros(3, N);
x_true(:,1) = [theta_int(1); y_meas(1); 0];  % başlangıç tork 0 varsay
LQ = chol(Q, 'lower');
for i = 2:N
    dt_i = t(i) - t(i-1);
    x_true(:,i) = f_transition(x_true(:,i-1), dt_i, m, g, L, I, b, k) + LQ*randn(3,1);
end

% --- EKF kurulumu ---
x0 = [theta_int(1); y_meas(1); 0];   % [theta; omega; tau]
ekf = extendedKalmanFilter( ...
    @(x) f_transition(x, t(2)-t(1), m, g, L, I, b, k), ...
    @h_measure, ...
    x0);

ekf.StateCovariance            = diag([10, 10, 10]);
ekf.ProcessNoise               = Q;
ekf.MeasurementNoise           = R;
ekf.StateTransitionJacobianFcn = @(x) F_jacobian(x, t(2)-t(1), m, g, L, I, b, k);
ekf.MeasurementJacobianFcn     = @H_jacobian;

% --- EKF döngüsü ---
x_est = zeros(3, N);
x_est(:,1) = x0;
for i = 2:N
    dt_i = t(i) - t(i-1);
    ekf.StateTransitionFcn         = @(x) f_transition(x, dt_i, m, g, L, I, b, k);
    ekf.StateTransitionJacobianFcn = @(x) F_jacobian(x, dt_i, m, g, L, I, b, k);
    predict(ekf);                    % model tahmini
    correct(ekf, y_meas(i));         % ölçüm düzeltmesi (omega)
    x_est(:,i) = ekf.State;
end

% --- Post-hoc: net tork doğrulama (kıyas amaçlı) ---
% EKF’den omega ve theta ile sayısal alpha:
alpha_num = [0, diff(x_est(2,:))./max(eps, diff(t)')];
tau_post  = I*alpha_num + m*g*L*sin(x_est(1,:)) + b*x_est(2,:) + k*x_est(1,:);

% --- Grafikler ---
figure('Position',[800 120 950 750]);

subplot(4,1,1);
plot(t, rad2deg(x_true(1,:)), 'b-', 'DisplayName','Gerçek θ (deg)'); hold on;
plot(t, rad2deg(x_est(1,:)) , 'r--','DisplayName','EKF θ (deg)');
plot(t, rad2deg(theta_int)  , 'g:','DisplayName','Gyro İntegral (deg)');
ylabel('Açı (deg)'); title('Açı (θ)'); legend; grid on;

subplot(4,1,2);
plot(t, rad2deg(x_true(2,:)), 'b-', 'DisplayName','Gerçek ω (deg/s)'); hold on;
plot(t, rad2deg(x_est(2,:)) , 'r--','DisplayName','EKF ω (deg/s)');
plot(t, rad2deg(y_meas)     , 'g:','DisplayName','Ölçüm ω (deg/s)');
ylabel('Açısal hız'); title('Açısal Hız (ω)'); legend; grid on;

subplot(4,1,3);
plot(t, rad2deg(alpha_kin), 'm-', 'DisplayName','α_kin = a_t/L (deg/s^2)'); hold on;
plot(t, rad2deg(alpha_num), 'k--','DisplayName','α_num (EKF türev) (deg/s^2)');
ylabel('Açısal ivme'); title('Açısal İvme Karşılaştırması'); legend; grid on;

subplot(4,1,4);
plot(t, x_est(3,:), 'r--', 'DisplayName','EKF τ (Nm)'); hold on;
plot(t, tau_post ,  'k:',  'DisplayName','Post-hoc τ (Nm)');
ylabel('Tork (Nm)'); xlabel('Zaman (s)'); title('Tork Tahmini'); legend; grid on;

% --- Basit hata metrikleri (simülasyon varsa anlamlı) ---
th_err = x_true(1,:) - x_est(1,:);
om_err = x_true(2,:) - x_est(2,:);
tau_err= x_true(3,:) - x_est(3,:);
fprintf('θ MSE   : %.4g rad^2\n', mean(th_err.^2));
fprintf('ω MSE   : %.4g (rad/s)^2\n', mean(om_err.^2));
fprintf('τ MSE   : %.4g Nm^2\n', mean(tau_err.^2));

%% ---------------- Yerel Fonksiyonlar (dosyanın en sonunda) ----------------
function xNext = f_transition(x, dt, m, g, L, I, b, k)
    % x = [theta; omega; tau]
    theta = x(1); omega = x(2); tau = x(3);
    omega_dot = (tau - m*g*L*sin(theta) - b*omega - k*theta) / I;
    theta_next = theta + omega*dt;
    omega_next = omega + omega_dot*dt;
    % tau'yu yavaş değişen süreç olarak modelle (random walk)
    tau_next = tau;  % gürültü Q(3,3) ile gelecek
    xNext = [theta_next; omega_next; tau_next];
end

function y = h_measure(x)
    % Ölçüm: sadece omega
    y = x(2);
end

function F = F_jacobian(x, dt, m, g, L, I, b, k)
    theta = x(1); % omega = x(2); tau = x(3);
    % df/dx
    d_omegadot_dtheta = ( -m*g*L*cos(theta) - k ) / I;
    d_omegadot_domega = ( -b ) / I;
    d_omegadot_dtau   = 1 / I;

    F = [ 1,              dt,                0;
          d_omegadot_dtheta*dt, 1 + d_omegadot_domega*dt, d_omegadot_dtau*dt;
          0,              0,                 1 ];
end

function H = H_jacobian(~)
    % y = [0 1 0] * x
    H = [0 1 0];
end

