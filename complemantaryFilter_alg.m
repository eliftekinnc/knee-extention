%The complementaryFilter parameters AccelerometerGain and MagnetometerGain
% can be tuned to change the amount each that the measurements of each 
% sensor impact the orientation estimate. When AccelerometerGain is set 
% to 0, only the gyroscope is used for the x- and y-axis orientation. 
% When AccelerometerGain is set to 1, only the accelerometer is used for 
% the x- and y-axis orientation. When MagnetometerGain is set to 0, only 
% the gyroscope is used for the z-axis orientation. When MagnetometerGain 
% is set to 1, only the magnetometer is used for the z-axis orientation.




clear; close all; clc;
load('knee_data.mat');
g     = 9.81;

a_n = knee.ACCYG(1:4967)*g;
a_t = knee.ACCZG(1:4967)*g;
a_x = knee.ACCXG(1:4967)*g;

w = knee.GYROXdegs(1:4967)*pi/180;
w_y = knee.GYROYdegs(1:4967)*pi/180;
w_z = knee.GYROZdegs(1:4967)*pi/180;

acc = [a_x a_n a_t];
gyro = [w w_y w_z];
pp = poseplot;

xlabel('x');
ylabel('y');
zlabel('z');

quat_log = quaternion.zeros(size(acc,1),1);

roll0 = 0;
pitch0 = 0;
yaw0 = 0;

% Disable magnetometer input.
cfilt = complementaryFilter(SampleRate=148.1481 ,HasMagnetometer=false);

for ii=1:size(acc,1)
    qimu = cfilt(acc(ii,:),gyro(ii,:))
    set(pp,"Orientation",qimu)
    drawnow limitrate
    %pause(0.01)
    quat_log(ii) = qimu;
end


euler_angles = eulerd(quat_log, "ZYX", "frame");

plot(euler_angles);
legend('Yaw','Pitch','Roll');
xlabel('Örnek');
ylabel('Açı (°)');
title('Complemantary Filtre Quaternion → Euler Dönüşümü');



