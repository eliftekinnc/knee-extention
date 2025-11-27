function [omega_hat, K_hist, P_hist] = kalman1d_gyro(z, Q, R, x0, P0)
% KALMAN1D_GYRO : En basit 1D Kalman filtresi (durum = omega)
% z  : gyro ölçümü (Nx1, rad/s)
% Q  : süreç gürültüsü (skaler). Örn: (alpha_std*dt)^2  (bkz. notlar)
% R  : ölçüm gürültüsü (skaler). Örn: var(z_static)
% x0 : başlangıç omega tahmini (varsayılan: z(1))
% P0 : başlangıç kovaryansı (varsayılan: 1)
% Çıktılar:
%   omega_hat : filtrelenmiş açısal hız (Nx1)
%   K_hist    : Kalman kazancı (Nx1)
%   P_hist    : hata kovaryansı (Nx1)

z = z(:);                % sütun vektör
N = numel(z);

if nargin < 4 || isempty(x0), x0 = z(1); end
if nargin < 5 || isempty(P0), P0 = 1;    end

omega_hat = zeros(N,1);
K_hist    = zeros(N,1);
P_hist    = zeros(N,1);

x = x0;         % durum: omega
P = P0;

for k = 1:N
    % --- Tahmin (predict) ---
    x_pred = x;           % F = 1
    P_pred = P + Q;       % P^- = P + Q

    % --- Güncelle (update) ---
    % H = 1  -> S = P^- + R,  K = P^-/(P^-+R)
    S = P_pred + R;
    K = P_pred / S;
    y = z(k) - x_pred;    % inovasyon
    x = x_pred + K*y;     % x = x^- + K*y
    P = (1 - K) * P_pred; % P = (I-KH)P^-

    % kayıt
    omega_hat(k) = x;
    K_hist(k)    = K;
    P_hist(k)    = P;
end
end