clear all; clc;
load("kinovea1.mat");
load("veri1.mat");
g     = 9.81;
fs = 148;


% IMU 1 (THIGH)

%5 numaralı sensor ÜST BACAK
a_x = veri1.ACCX_G_2*g;
a_y = veri1.ACCY_G_2*g;
a_z = veri1.ACCZ_G_2*g;

w_x = veri1.GYROX_deg_s_2*pi/180;
w_y = veri1.GYROY_deg_s_2*pi/180;
w_z = veri1.GYROZ_deg_s_2*pi/180;

acc = [a_x a_y a_z];
gyro = [w_x w_y w_z];








% IMU 2 (SHANK)

% %2 numaralı sensor BACAĞIN YANI VE OK İLERİ (pitch)
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
% 
% 


% 4 nuamralı sensor YAN VE OK YUKARI (pitch)
% a_x1 = veri1.ACCX_G_*g;
% a_y1 = veri1.ACCY_G_*g;
% a_z1 = veri1.ACCZ_G_*g;
% 
% w_x1 = veri1.GYROX_deg_s_*pi/180;
% w_y1 = veri1.GYROY_deg_s_*pi/180;
% w_z1 = veri1.GYROZ_deg_s_*pi/180;
% 




% 8 numaralı sensor ARKA BACAK VE OK AŞAĞI
a_x1 = veri1.ACCX_G_3*g;
a_y1 = veri1.ACCY_G_3*g;
a_z1 = veri1.ACCZ_G_3*g;

w_x1 = veri1.GYROX_deg_s_3*pi/180;
w_y1 = veri1.GYROY_deg_s_3*pi/180;
w_z1 = veri1.GYROZ_deg_s_3*pi/180;



acc1 = [a_x1 a_y1 a_z1];
gyro1 = [w_x1 w_y1 w_z1];

N = size(acc1,1);


[q_results, angles] = mahony_filter(acc, gyro, fs); %üst bacak quat

[q_results1, angles1] = mahony_filter(acc1, gyro1, fs); %alt bacak quat

q_thigh = quaternion(q_results);
q_shank = quaternion(q_results1);

q_rel = quaternion.zeros(N,1);

for a=1:1332
    q_thigh0 = q_thigh(a);
    q_shank0 = q_shank(a);
    
    for i=1:N
        
        q_thigh_corr(i) = conj(q_thigh0) * q_thigh(i);
        q_shank_corr(i) = conj(q_shank0) * q_shank(i);
    
    end
    
    q_rel = conj(q_thigh_corr) .* q_shank_corr;
    
    [w,x,y,z] = parts(q_rel);
    
    w = max(min(w,1),-1);
    theta = 2*acos(w);
    axis = zeros(N,3);
    
    for i=1:N
        
        s = sin(theta(i)/2);
        
        if abs(s) > 1e-6
            axis(i,:) = [x(i) y(i) z(i)]/s;
        end
        
    end
    
    [coeff,~,~] = pca(axis);
    
    knee_axis = coeff(:,1);
    knee_axis = knee_axis / norm(knee_axis);
    
    knee_angle = zeros(N,1);
    
    for i=1:N
        
        knee_angle(i) = theta(i) * (axis(i,:) * knee_axis);
        
    end
    
    knee_angle = -rad2deg(knee_angle);
    
    knee_angle = detrend(knee_angle);
    
    
    
    knee_angle(abs(knee_angle) > 120) = NaN;
    knee_angle = fillmissing(knee_angle,'linear');
    % figure;
    % plot(knee_angle+22);
    
    
    % figure;
    % plot(kinovea_data.A__1-113);
    % 
    
    
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
    
    rmse = sqrt(mean(hata3.^2))
    a
end 
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

