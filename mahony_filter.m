function [q, euler] = mahony_filter(accel, gyro, fs)
    % accel: [N x 3] ivmeölçer verisi (m/s^2)
    % gyro:  [N x 3] jiroskop verisi (rad/s)
    % fs:    Örnekleme frekansı (Hz)

    % --- Parametreler (C kodundaki tanımlar) ---
    Kp = 2.0;       % Oransal kazanç
    Ki = 0.005;     % İntegral kazanç
    halfT = 1/(2*fs); % Örnekleme periyodunun yarısı

    N = size(accel, 1);
    q = zeros(N, 4);
    euler = zeros(N, 3);

    % Başlangıç değerleri
    q0 = 1; q1 = 0; q2 = 0; q3 = 0;
    exInt = 0; eyInt = 0; ezInt = 0;

    for i = 1:N
        gx = gyro(i,1); gy = gyro(i,2); gz = gyro(i,3);
        ax = accel(i,1); ay = accel(i,2); az = accel(i,3);

        % 1. İvmeölçer normalizasyonu
        norm = sqrt(ax*ax + ay*ay + az*az);
        if norm > 0
            ax = ax / norm; ay = ay / norm; az = az / norm;
        end

        % 2. Tahmini yerçekimi yönü (v)
        vx = 2*(q1*q3 - q0*q2);
        vy = 2*(q0*q1 + q2*q3);
        vz = q0*q0 - q1*q1 - q2*q2 + q3*q3;

        % 3. Hata hesaplama (Cross product: a x v)
        ex = (ay*vz - az*vy);
        ey = (az*vx - ax*vz);
        ez = (ax*vy - ay*vx);

        % 4. İntegral hata birikimi
        exInt = exInt + ex * Ki;
        eyInt = eyInt + ey * Ki;
        ezInt = ezInt + ez * Ki;

        % 5. Jiroskop verisini düzeltme
        gx = gx + Kp*ex + exInt;
        gy = gy + Kp*ey + eyInt;
        gz = gz + Kp*ez + ezInt;

        % 6. Kuaterniyon güncelleme (Entegrasyon)
        q0_new = q0 + (-q1*gx - q2*gy - q3*gz)*halfT;
        q1_new = q1 + (q0*gx + q2*gz - q3*gy)*halfT;
        q2_new = q2 + (q0*gy - q1*gz + q3*gx)*halfT;
        q3_new = q3 + (q0*gz + q1*gy - q2*gx)*halfT;
        
        % Değerleri güncelle
        q0 = q0_new; q1 = q1_new; q2 = q2_new; q3 = q3_new;

        % 7. Kuaterniyon normalizasyonu
        norm_q = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
        q0 = q0 / norm_q; q1 = q1 / norm_q; q2 = q2 / norm_q; q3 = q3 / norm_q;

        % Sonuçları kaydet
        q(i, :) = [q0, q1, q2, q3];
        
        % ZYX sırasıyla Euler açılarına çevir (derece)
        euler(i, :) = quat2eul([q0, q1, q2, q3], 'ZYX') * 180/pi;
    end
end