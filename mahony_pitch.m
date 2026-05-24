clear all; clc;
load("kinovea1.mat");
load("veri1.mat");
g     = 9.81;


% 5 numaralı sensor ÜST BACAK
% a_x1 = veri1.ACCX_G_2*g;
% a_y1 = veri1.ACCY_G_2*g;
% a_z1 = veri1.ACCZ_G_2*g;
% 
% w_x1 = veri1.GYROX_deg_s_2*pi/180;
% w_y1 = veri1.GYROY_deg_s_2*pi/180;
% w_z1 = veri1.GYROZ_deg_s_2*pi/180;





% 2 numaralı sensor BACAĞIN YANI VE OK İLERİ (pitch)
% a_x1 = veri1.ACCX_G_5*g;
% a_y1 = veri1.ACCY_G_5*g;
% a_z1 = veri1.ACCZ_G_5*g;
% 
% w_x1 = veri1.GYROX_deg_s_5*pi/180;
% w_y1 = veri1.GYROY_deg_s_5*pi/180;
% w_z1 = veri1.GYROZ_deg_s_5*pi/180;
% 



%7 numaralı sensor BACAĞIN ÖNÜ VE OK ÇAPRAZ (roll)
% a_x1 = veri1.ACCX_G_1*g;
% a_y1 = veri1.ACCY_G_1*g;
% a_z1 = veri1.ACCZ_G_1*g;
% 
% w_x1 = veri1.GYROX_deg_s_1*pi/180;
% w_y1 = veri1.GYROY_deg_s_1*pi/180;
% w_z1 = veri1.GYROZ_deg_s_1*pi/180;
% 




% 4 nuamralı sensor YAN VE OK YUKARI (pitch)
a_x1 = veri1.ACCX_G_*g;
a_y1 = veri1.ACCY_G_*g;
a_z1 = veri1.ACCZ_G_*g;

w_x1 = veri1.GYROX_deg_s_*pi/180;
w_y1 = veri1.GYROY_deg_s_*pi/180;
w_z1 = veri1.GYROZ_deg_s_*pi/180;





% 8 numaralı sensor ARKA BACAK VE OK AŞAĞI
% a_x1 = veri1.ACCX_G_3*g;
% a_y1 = veri1.ACCY_G_3*g;
% a_z1 = veri1.ACCZ_G_3*g;
% 
% w_x1 = veri1.GYROX_deg_s_3*pi/180;
% w_y1 = veri1.GYROY_deg_s_3*pi/180;
% w_z1 = veri1.GYROZ_deg_s_3*pi/180;



acc = [a_x1 a_y1 a_z1];
gyro = [w_x1 w_y1 w_z1];

fs = 148; % Örnekleme hızı

[q_results, angles] = mahony_filter(acc, gyro, fs);

% figure;
% plot(angles(:, 2)-2); % 2. sütun Pitch açısı
% title('Diz Bükülme Açısı (Flexion/Extension)');
% ylabel('Derece');
% xlabel('Örnek Sayısı');
% 
% 
% figure;
% plot(kinovea_data.A__1-113);
% 

q_results = quaternion(q_results);
euler_angles1 = eulerd(q_results, "ZYX", "frame");
knee_angle = euler_angles1(:,2);
knee_angle  = rad2deg(unwrap(deg2rad(knee_angle)));


% plot(euler_angles1);
% legend('Yaw','Pitch','Roll');
% xlabel('Örnek');
% ylabel('Açı (°)');
% title('IMU Filtre Quaternion → Euler Dönüşümü');

% figure;
% subplot(2,1,1); plot(imu_data.ACCXTimeSeries_s_, euler_angles1);
% subplot(2,1,2); plot(kinovea_data.Zaman_ms_/1000, kinovea_data.A__1);
% 
% figure;
% subplot(2,1,1); plot(imu_data.ACCXTimeSeries_s_(450:end), euler_angles1(450:end, 2));
% subplot(2,1,2); plot(kinovea_data.Zaman_ms_(243:end)/1000, kinovea_data.A__1(243:end));





%  Zaman dizileri
t_imu = veri1.ACCXTimeSeries_s_(1:end);
t_kinovea = (kinovea_angle_1.Zaman_ms_)/ 1000;

% Ham veriler
val_method3 = knee_angle-mean(knee_angle(1:148*10));
val_kinovea = kinovea_angle_1.A__1-mean(kinovea_angle_1.A__1(1:30*10));

% INTERPOLATION
kinovea_resampled = interp1(t_kinovea, val_kinovea, t_imu, 'pchip');

[xa, ya, D] = alignsignals(val_method3, kinovea_resampled);

% Ortak uzunluk
L = min([length(t_imu) length(xa) length(ya)]);


t = t_imu(D:L)-t_imu(D);
imu = xa(D:L);
kino = ya(D:L);

hata3 = imu - kino;

rmse = sqrt(mean(hata3.^2));

% 1. Ana Figür
figure('Color','w','Name','MatlabEKF+Euler vs Kinovea Karşılaştırması');


subplot(2,1,1)

plot(t, imu,'r','LineWidth',1.5)
hold on
plot(t, kino,'b','LineWidth',1.5)

title(['Diz Açısı Kıyaslaması (RMSE: ', num2str(rmse,'%.2f'),'°)'])
ylabel('Açı (Derece)')
legend('MatlabEKF+Euler','Kinovea')
grid on


subplot(2,1,2)

fill([t fliplr(t)], [hata3 fliplr(zeros(size(hata3)))], ...
     [0.9 0.9 0.9], 'EdgeColor','none')
hold on

plot(t,hata3,'k','LineWidth',1.2)

title('Anlık Hata Dağılımı')
xlabel('Zaman (saniye)')
ylabel('Fark (°)')
