load("kinovea1.mat");
load("veri1.mat");
g=9.81;
fs = 148;



% 5 numaralı sensor ÜST BACAK
% a_x1 = veri1.ACCX_G_2*g;
% a_y1 = veri1.ACCY_G_2*g;
% a_z1 = veri1.ACCZ_G_2*g;
% 
% w_x1 = veri1.GYROX_deg_s_2*pi/180;
% w_y1 = veri1.GYROY_deg_s_2*pi/180;
% w_z1 = veri1.GYROZ_deg_s_2*pi/180;





% 2 numaralı sensor BACAĞIN YANI VE OK İLERİ
% a_x1 = veri1.ACCX_G_5*g;
% a_y1 = veri1.ACCY_G_5*g;
% a_z1 = veri1.ACCZ_G_5*g;
% 
% w_x1 = veri1.GYROX_deg_s_5*pi/180;
% w_y1 = veri1.GYROY_deg_s_5*pi/180;
% w_z1 = veri1.GYROZ_deg_s_5*pi/180;




%7 numaralı sensor BACAĞIN ÖNÜ VE OK ÇAPRAZ
% a_x1 = veri1.ACCX_G_1*g;
% a_y1 = veri1.ACCY_G_1*g;
% a_z1 = veri1.ACCZ_G_1*g;
% 
% w_x1 = veri1.GYROX_deg_s_1*pi/180;
% w_y1 = veri1.GYROY_deg_s_1*pi/180;
% w_z1 = veri1.GYROZ_deg_s_1*pi/180;





% 4 nuamralı sensor YAN VE OK YUKARI
% a_x1 = veri1.ACCX_G_*g;
% a_y1 = veri1.ACCY_G_*g;
% a_z1 = veri1.ACCZ_G_*g;
% 
% w_x1 = veri1.GYROX_deg_s_*pi/180;
% w_y1 = veri1.GYROY_deg_s_*pi/180;
% w_z1 = veri1.GYROZ_deg_s_*pi/180;
% 




%8 numaralı sensor ARKA BACAK VE OK AŞAĞI
a_x1 = veri1.ACCX_G_3*g;
a_y1 = veri1.ACCY_G_3*g;
a_z1 = veri1.ACCZ_G_3*g;

w_x1 = veri1.GYROX_deg_s_3*pi/180;
w_y1 = veri1.GYROY_deg_s_3*pi/180;
w_z1 = veri1.GYROZ_deg_s_3*pi/180;


acc = [a_x1 a_y1 a_z1];
gyro = [w_x1 w_y1 w_z1];

quat_log = quaternion.zeros(size(acc,1),1);

ifilt = imufilter(SampleRate=148);
for ii=1:size(acc,1)
    qimu = ifilt(acc(ii,:),gyro(ii,:));
    %set(pp,"Orientation",qimu)
    %drawnow limitrate
    quat_log(ii) = qimu;
end

quat_log_conj = conj(quat_log);

N = size(acc,1);

%delta quaternions
dq = quaternion.zeros(size(acc,1)-1,1);

for i=1:N-1
    dq(i) = quat_log(i+1)*quat_log_conj(i);
end

% axis hesaplama

[w, x, y, z] = parts(dq);


N = length(x); 

theta = 2*acos(w);

axis = [];
for i = 1:N

    s = sin(theta(i)/2);

    if abs(s) > 1e-6
        axis= [axis; x(i)/s y(i)/s z(i)/s] ;
    end

end





[coeff, score, latent] = pca(axis);

ana_eksen = coeff(:, 1);


% ana_eksen'in [1x3] satır vektörü olduğundan emin olalım
if iscolumn(ana_eksen)
    ana_eksen = ana_eksen'; % [1x3] satır vektörüne çevir
end


ham_eksenler = [x, y, z];

dot_products = ham_eksenler * ana_eksen';

yon_isaretleri = sign(dot_products);

% 4. Gerçek (Yönlü) Açıları Hesaplayın
% Büyüklük * İşaret
anlik_acilar = theta .* yon_isaretleri; % Radyan cinsinden

anlik_acilar_derece = rad2deg(anlik_acilar);

anlik_acilar_derece(abs(anlik_acilar_derece) > 5) = 0;

% figure;
% plot(cumtrapz(anlik_acilar_derece));
% title('dizin ekstansiyon-fleksiyon açısı');
% xlabel('zaman');
% ylabel('Açı');
% grid on;

ana_eksen_n = ana_eksen/norm(ana_eksen);

N_sample = size(anlik_acilar_derece, 1);

fi = zeros(N_sample, 1);
%anlık açıların ana eksen üzerindeki iz düşümü 
for i=1:N_sample
    ham_eksenler(i,:) = ham_eksenler(i,:)/norm(ham_eksenler(i,:));
    fi(i) = (ham_eksenler(i,:)*ana_eksen_n') .* abs(anlik_acilar_derece(i));
end


knee_angle = -cumtrapz(fi);
% figure;
% plot(knee_angle);
% title('dizin ekstansiyon-fleksiyon açısının ana eksene izdüşümü');
% xlabel('zaman');
% ylabel('Açı');
% grid on;



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
figure('Color','w','Name','MatlabEKF+PCA vs Kinovea Karşılaştırması');


subplot(2,1,1)

plot(t, imu,'r','LineWidth',1.5)
hold on
plot(t, kino,'b','LineWidth',1.5)

title(['Diz Açısı Kıyaslaması (RMSE: ', num2str(rmse,'%.2f'),'°)'])
ylabel('Açı (Derece)')
legend('MatlabEKF+PCA','Kinovea')
grid on


subplot(2,1,2)

fill([t fliplr(t)], [hata3 fliplr(zeros(size(hata3)))], ...
     [0.9 0.9 0.9], 'EdgeColor','none')
hold on

plot(t,hata3,'k','LineWidth',1.2)

title('Anlık Hata Dağılımı')
xlabel('Zaman (saniye)')
ylabel('Fark (°)')
grid on