fs = 1259,2593;
t = knee.EMG1TimeSeriess;
tseri = knee.ACCZTimeSeriess(1:4967);
% EMG verileri
vastus = knee.EMG1mV(:); % Vastus (quadriceps)
rectus = knee.EMG1mV1(:); % Rectus femoris
gastro = knee.EMG1mV2(:); % Gastrocnemius


% Filtreleme ve zarf tespiti
vastus_f = medfilt1(vastus, 15);
rectus_f = medfilt1(rectus, 15);
gastro_f = medfilt1(gastro, 15);
vastus_env = lowpass(vastus_f, 10, fs);
rectus_env = lowpass(rectus_f, 10, fs);
gastro_env = lowpass(gastro_f, 10, fs);

% Normalizasyon
vastus_norm = vastus_env / max(vastus_env);
rectus_norm = rectus_env / max(rectus_env);
gastro_norm = gastro_env / max(gastro_env);

% Kuvvet tahmini
k_vastus = 50; % Kuvvet katsayısı (deneysel)
k_rectus = 50;
k_gastro = 30;
F_vastus = k_vastus * vastus_norm;
F_rectus = k_rectus * rectus_norm;
F_gastro = k_gastro * gastro_norm;

% Tork tahmini
L_gastro = 0.25; % Gastrocnemius moment kolu
tau_emg = (F_vastus * L * 1) + (F_rectus * L * 1) - (F_gastro * L_gastro * 1);

% Füzyon
w_emg = 0.7; % EMG ağırlığı
w_post = 0.3; % Post-hoc ağırlığı
tau_fusion = w_emg * tau_emg + w_post * tau_post;

% Mevcut grafik güncellemesi
subplot(3,1,3);
yyaxis left
plot(tseri, rad2deg(alpha_num), 'm-', 'DisplayName', '\alpha_{num} (deg/s^2)');
yyaxis right
plot(t, tau_fusion, 'g-', 'DisplayName', '\tau_{fusion} (Nm)');
hold on;
plot(t, tau_post, 'k-', 'DisplayName', '\tau_{post} (Nm)');
plot(t, tau_emg, 'b--', 'DisplayName', '\tau_{EMG} (Nm)');
xlabel('Zaman (s)'); title('Açısal İvme ve Tork'); legend; grid on;