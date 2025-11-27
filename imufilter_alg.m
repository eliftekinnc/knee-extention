%The imufilter parameters are a subset of the ahrsfilter parameters. 
% The AccelerometerNoise, GyroscopeNoise, MagnetometerNoise, and 
% GyroscopeDriftNoise are measurement noises. The sensors' datasheets help
% determine those values.


clear; close all; clc;
load('knee_data.mat');
g     = 9.81;


a_n = knee.ACCYG(1:4967)*g;
a_t = knee.ACCZG(1:4967)*g;
a_x = knee.ACCXG(1:4967)*g;

w = knee.GYROXdegs(1:4967)*pi/180;
w_y = knee.GYROYdegs(1:4967)*pi/180;
w_z = knee.GYROZdegs(1:4967)*pi/180;

sifir = zeros(4967,1);

acc = [a_x a_n a_t];
gyro = [w w_y w_z];
pp = poseplot;

xlabel('x');
ylabel('y');
zlabel('z');

quat_log = quaternion.zeros(size(acc,1),1);


% Disable magnetometer input.
ifilt = imufilter(SampleRate=148);
for ii=1:size(acc,1)
    qimu = ifilt(acc(ii,:),gyro(ii,:))
    set(pp,"Orientation",qimu)
    drawnow limitrate
    %pause(0.01)

    quat_log(ii) = qimu;
end

%quaternionları euler açılarına çevirelim

euler_angles = eulerd(quat_log, "ZYX", "frame");

plot(euler_angles);
legend('Yaw','Pitch','Roll');
xlabel('Örnek');
ylabel('Açı (°)');
title('IMU Filtre Quaternion → Euler Dönüşümü');

