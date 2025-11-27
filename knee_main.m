%% Diz EKF (extendedKalmanFilter obje ile, sim yok)
clear; close all; clc;
load('knee_data.mat');

% ----------------- Fizik parametreleri -----------------
g     = 9.81;
L     = 0.30;           
L_tam = 0.50;            
m     = 5.00;           
I     = (1/3)*m*L_tam^2; 

% ----------------- Veri boyutu -----------------
N = min([4967, numel(knee.ACCXTimeSeriess), numel(knee.GYROXdegs), numel(knee.ACCYG), numel(knee.ACCZG)]);
t = knee.ACCXTimeSeriess(1:N);

% ----------------- Ham sinyaller -----------------
w_deg = (-1)*knee.GYROXdegs(1:N);   
a_t   = knee.ACCZG(1:N)*g;          

% ----------------- Filtre + bias -----------------
 w_f = medfilt1(w_deg,15);
 a_t_f = medfilt1(a_t,15);

K0 = min(150, floor(N/10));
w_f   = w_f   - mean(w_f(1:K0));
a_t_f = a_t_f - mean(a_t_f(1:K0));

% ----------------- Radyan + referans -----------------
w_rad   = deg2rad(w_f);        
theta_i = cumtrapz(t, w_rad);   

% ----------------- Ölçüm (gyro) -----------------
y = w_rad;                      

% ----------------- EKF (object) -----------------
% Durum: x = [theta; omega]
Q = diag([1e-3, 1e-2]);         % süreç gürültüsü
R = (0.05)^2;                   % ölçüm gürültüsü var(y) ≈ (rad/s)^2
x0 = [theta_i(1); y(1)];

% Başlat (dt geçici; her döngüde güncellenecek)
ekf = extendedKalmanFilter( ...
    @(x) f_step(x, t(2)-t(1), m, g, L, I), ...
    @h_meas, x0);

ekf.StateCovariance            = eye(2)*10;
ekf.ProcessNoise               = Q;
ekf.MeasurementNoise           = R;
ekf.StateTransitionJacobianFcn = @(x) F_jac(x, t(2)-t(1), m, g, L, I);
ekf.MeasurementJacobianFcn     = @H_jac;

x_est = zeros(2,N); x_est(:,1) = x0;

% inovasyon/NIS izleme
innov = zeros(N,1); Sarr = zeros(N,1);

for k = 2:N
    dt = t(k) - t(k-1);

    % adım-adım dt ile fonksiyonları güncelle
    ekf.StateTransitionFcn         = @(x) f_step(x, dt, m, g, L, I);
    ekf.StateTransitionJacobianFcn = @(x) F_jac(x, dt, m, g, L, I);

    % --- PREDICT ---
    predict(ekf);

    % inovasyon ve S (correct öncesi, predicted kovaryans ile)
    x_pred = ekf.State;
    P_pred = ekf.StateCovariance;
    z_pred = h_meas(x_pred);
    H = [0 1];
    S = H*P_pred*H' + R;

    innov(k) = y(k) - z_pred;
    Sarr(k)  = S;

    % --- UPDATE ---
    correct(ekf, y(k));

    x_est(:,k) = ekf.State;
end

% ----------------- Post-hoc tork (kontrol amaçlı) -----------------

omega_s  = movmean(x_est(2,:), 11);

alpha_num = [0, diff(omega_s)./max(eps, diff(t)')];     % rad/s^2
b_vis = 0.05;  k_el = 2.0;
tau_post = I*alpha_num + m*g*L*sin(x_est(1,:)) + b_vis*x_est(2,:) + k_el*x_est(1,:);

% NIS (ideal ortalama ≈ 1)
nis_mean = mean((innov.^2)./Sarr,'omitnan');
fprintf('NIS ~ %.2f (ideal ≈ 1)\n', nis_mean);

% ----------------- Grafikler -----------------
figure('Position',[900 100 900 720]);

subplot(3,1,1);
plot(t, rad2deg(x_est(1,:)),'r--','DisplayName','EKF \theta (deg)'); hold on;
plot(t, rad2deg(theta_i)   ,'b-','DisplayName','Gyro İntegral (deg)');
ylabel('Açı (deg)'); title('Diz Açısı'); legend; grid on;

subplot(3,1,2);
plot(t, rad2deg(x_est(2,:)),'m--','DisplayName','EKF \omega (deg/s)'); hold on;
plot(t, rad2deg(y)         ,'b:','DisplayName','Ölçülen gyro \omega (deg/s)');
ylabel('\omega (deg/s)'); title('Açısal Hız'); legend; grid on;

subplot(3,1,3);
yyaxis left
plot(t, rad2deg(alpha_num),'m-','DisplayName','\alpha_{num} (deg/s^2)'); ylabel('\alpha (deg/s^2)');
yyaxis right
plot(t, tau_post,'k-','DisplayName','\tau_{post} (Nm)'); ylabel('\tau (Nm)');
xlabel('Zaman (s)'); title('Açısal İvme ve Post-hoc Tork'); legend; grid on;

% ----------------- Yerel fonksiyonlar (dosyanın SONU) -----------------
function xNext = f_step(x, dt, m, g, L, I)
    % x = [theta; omega], u=0; omega_dot = -(m g L / I) sin(theta)
    th = x(1); om = x(2);
    om_dot = -(m*g*L/I)*sin(th);
    xNext = [th + om*dt;
             om + om_dot*dt];
end

function y = h_meas(x)
    % Ölçüm: y = omega
    y = x(2);
end

function F = F_jac(x, dt, m, g, L, I)
    th = x(1);
    F = [1, dt;
         -(m*g*L/I)*cos(th)*dt, 1];
end

function H = H_jac(~)
    H = [0 1];
end
