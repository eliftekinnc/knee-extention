g  = 9.81;
l  = 0.5;
l_a= 0.3;
m  = 5;
%% sampling frekansları farklı olduğu için istediğimiz verileri çekelim

w = (-1)*knee.GYROXdegs(1:4967)*pi/180;
a_n = (-1)*knee.ACCYG(1:4967)*g;
a_t = knee.ACCZG(1:4967)*g;
%% 
% *herhangi bir time series(emg hariç)*

tseri = knee.ACCXTimeSeriess(1:4967);
tstep = 0.00675;
%% 
% *BİASLERİ ÇIKARALIM*

w = w - mean(w(1:4967));
a_n = a_n - mean(a_n(1:4967));
a_t = a_t - mean(a_t(1:4967));
%% 
% *TETA ve ALFAYI TANIMLAYALIM*

teta = acos((w.^2*l_a - a_n)/g);
alfa = (a_t + g*sin(teta))/l_a;
%% 
% *KOMPLESK SAYILAR ÇIKTIĞI İÇİN*

abs(teta);
abs(alfa);
%% 
% *KONTROL EDELİM*

plot(tseri, abs(teta)) 
plot(tseri, abs(alfa)) 
%% 
% *DEĞERLER TUTARLI MI KONTROL EDELİM*

A_n = w.^2*l_a - g*cos(teta); %formulde yerine yazıyoruz
A_t = alfa*l_a - g*sin(teta);

error_1 = abs(mean((A_n - a_n).^2));
error_2 = abs(mean((A_t - a_t).^2));
%% 
% *alfa tetanın ikinci türevine karşılık geliyor mu? ne kadar geliyor?*

for i = 2:1:4966
    teta_turev1(i) = teta(i+1)-teta(i-1) / (2*tstep);
end

for i = 2:1:4965
    teta_turev2(i) = teta_turev1(i+1)-teta_turev1(i-1) / (2*tstep);
end

teta_turev2.';

mse = abs(mean( (teta_turev2.' - alfa(1:4965)).^2 )); % error çok büyük

%% 
% *MATEMATİKSEL MODELİ KULLANARAK BİR ŞEYLER DENEYELİM*

tork_o = (m * l^2/3) .* alfa + (m*g*l/2) .* sin(teta);

tseri2 = knee.EMG1TimeSeriess;
plot(tseri2, knee.EMG1mV);
plot(tseri2, knee.EMG1mV1);
plot(tseri2, knee.EMG1mV2);