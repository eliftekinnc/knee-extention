
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



p0 = [-1.1722, -9.6986, 0.9233]; %body frame
p1 = [0, 0, -9.8127]; %world frame

%normalize et
p0 = p0 / norm(p0);
p1 = p1 / norm(p1);


% cross ve dot product
c = cross(p0, p1) ; 
d = dot(p0, p1) ;


if d >= 1 - eps 
    R = eye(3);
elseif d<= -1 + eps

else
    teta = acos(d);
    a = [0 -c(3) c(2); c(3) 0 -c(1); -c(2) c(1) 0];
    R = eye(3) + sin(teta)*a + (1-cos(teta))*a*a;
end


%ivme verilerini hizalayalım
B = zeros(size(acc));  % sonucu saklamak için aynı boyutta matris

for i = 1:size(acc,1)
    row = acc(i,:);            % i. satırı al
    row_processed = R*row';  % örnek işlem: normalize et
    B(i,:) = row_processed;  % sonucu sakla
end

