PPGo = PPG3;                  % Select PPG data from one of the available subjects (e.g., PPG1 to PPG7)
PPG = PPGo(:,1);              % Reference PPG signal
PPGcont = PPGo(:,2);          % Contaminated PPG signal (affected by motion artifacts)
x = PPGo(:,3);                % Acceleration signal from X-axis
y = PPGo(:,4);                % Acceleration signal from Y-axis
z = PPGo(:,5);                % Acceleration signal from Z-axis

% Remove DC offset from all signals
PPG = PPG - mean(PPG);
PPGcont = PPGcont - mean(PPGcont);
x = x - mean(x);
y = y - mean(y);
z = z - mean(z);
n = sqrt( (x.^2) + (y.^2) + (z.^2) ); % Acceleration magnitude from the three axes

%%%%%%%%%%%%%%%%%%% Initial parameters for adaptive filter %%%%%%%%%%%%%%%%%%%%%%%
k = length(PPG); 
N = 200;                 % Filter order
M = 10;                  % Projection order

miu = 0.03;              % Step size (mu), should be less than 2 for convergence
gamma = 0.002;           % Regularization parameter to prevent numerical instability
C = 0.002;               % Positive constant used in variable step-size algorithms
alpha = 0.2;             % Forgetting factor (0 < alpha < 1), used in VSSAPA

mov = [x y z n y z]; % Matrix of motion signals: x, y, z, n, xy, xyz
R = zeros(k,6);      % Matrix to store recovered signals: [Rx; Ry; Rz; Rn; Rxy; Rxyz]

for j=1:6
    W=zeros(N,1);        % Weight vector (filter coefficients)
    e=zeros(1,k);        % Error signal
    X=zeros(N,M);        % Matrix for projections (input data for the filter)
    P=zeros(N,1);        % Placeholder vector (used in VSSNLMS)
    r = mov(:,j);        % Select the motion signal to be used as the reference input
    
    % Select the appropriate noisy signal based on the current motion component
    if (j == 5)
        noisy = R(:, 1);   % For Rxy, use the result from filtering with x-axis (Rx)
    elseif (j == 6)
        noisy = R(:, 5);   % For Rxyz, use the result from filtering with Rxy
    else
        noisy = PPGcont;   % For Rx, Ry, Rz, and Rn use the original contaminated PPG
    end
    
% Adaptive filtering loop
    for i = N:k
        delay = r(i:-1:i-N+1);           % Delayed input vector
        X = [delay X(:,1:M-1)];          % Matrix of input projections
        Y = X'*W;                        % Output of the adaptive filter
        E = noisy(i:-1:i-M+1)-Y;         % Error signal (difference between actual noisy signal and filter output)
        e(i) = mean(E);                  % Average error for this sample

        % === Select ONE adaptive algorithm (uncomment only one) ===
        W = apa(W,miu,X,E,M,gamma);                % APA
        % W = vapa(W,miu,X,E);                     % VAPA
        % W = vssapa(W,miu,X,E,C,alpha,P,M,gamma); % VSSAPA
    end
        R(:,j) = e';     % Error signal (filtered output)
end

t = ( 0 : length( PPG ) - 1 ) / 100 ;  % Discrete time vector

figure(1)
plot(t, PPG, 'LineWidth', 1.7); hold on
plot(t, PPGcont + 6, 'LineWidth', 1.7); hold on

% Assuming R has the same number of rows as PPG
plot(t, R(:,1) + 1, 'LineWidth', 1.7); hold on
plot(t, R(:,2) + 2, 'LineWidth', 1.7); hold on
plot(t, R(:,3) + 3, 'LineWidth', 1.7); hold on
plot(t, R(:,4) + 4, 'LineWidth', 1.7); hold on
plot(t, R(:,6) + 5, 'LineWidth', 1.7); grid on

legend({'PPG ref', 'PPG cont', 'Rx', 'Ry', 'Rz', 'Rn', 'Rxyz'})
xlabel('Tiempo (s)'); ylabel('Amplitud (v)');
% xlim([82.13 115.92]); %ylim([-0.5 7])

R = R(:,2);         % Best recovery (from filtered signals)
ruido = y ;         % Motion signal used to filter PPG signal
RR = R - mean(R);   % Remove DC offset from filtered signal

%%%%%%%%%%%%%%%%%%%%%%%%% Power Spectrum (PSD) %%%%%%%%%%%%%%%%%%%%%%%%%%%
window = boxcar(128);     % Rectangular window
noverlap = 64;            % 50% overlap
nfft = 512;               % Section size
fs = 100;                 % Sampling frequency in Hz

[PSD_PPG, f_PPG] = pwelch(PPG, window, noverlap, nfft, fs); % PSD of reference PPG signal
[PSD_PPGcont, f_PPGcont] = pwelch(PPGcont, window, noverlap, nfft, fs); % PSD of contaminated PPG signal
[PSD_R, f_R] = pwelch(RR, window, noverlap, nfft, fs); % PSD of filtered signal
[PSD_r, f_r] = pwelch(ruido, window, noverlap, nfft, fs); % PSD of noise

PSD_PPGA = PSD_PPG * (max(PSD_PPGcont(4:end))/max(PSD_PPG(4:end))) ;

figure(2)
ax1 = axes('Position',[0.2 0.1 0.7 0.8]);
plot(ax1, f_PPG,PSD_PPGA,f_PPGcont,PSD_PPGcont,f_R,PSD_R,'LineWidth',2); 
xlim(ax1, [0 10]); % ylim(ax1, [0 0.0045]);
legend(ax1, 'Pulso de referencia', 'Pulso contaminado','Pulso recuperado')
xlabel(ax1, 'Frecuencia (Hz)'); ylabel(ax1, 'PSD (mV^2/Hz)'); grid on
title(ax1, 'Densidad espectral de potencia')

ax2 = axes('Position',[0.52 0.28 0.35 0.5]); % [0.64 0.28 0.25 0.5] ); %
plot(ax2, f_PPG,PSD_PPGA,'LineWidth',2); hold on; 
plot(ax2, f_PPGcont,PSD_PPGcont,'LineWidth',2); hold on
plot(ax2, f_R,PSD_R,'LineWidth',2); grid on
xlim(ax2, [2.34 7]); % ylim(ax2, [0 0.006]);
xlabel(ax2, 'Frecuencia (Hz)'); ylabel(ax2, 'PSD (mV^2/Hz)'); 

figure(3)
plot(f_r,PSD_r,'LineWidth',1.7,'Color',[32/225,178/225,170/225]); 
xlim([0 10]); grid on; legend('Señal del acelerómetro')
xlabel('Frecuencia (Hz)'); ylabel('PSD (mV^2/Hz)'); 
title('Densidad espectral de potencia')

%%%%%%%%%%%%%%%%%%%%%% Spectral Coherence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[C_P_R, F_P_R] = mscohere(PPG, RR, window, noverlap, nfft, fs);         % Coherence: PPG vs. recovered
[C_Pcont_R, F_Pcont_R] = mscohere(PPGcont, RR, window, noverlap, nfft, fs); % Coherence: PPGcont vs. recovered
[C_P_Pcont, F_P_Pcont] = mscohere(PPG, PPGcont, window, noverlap, nfft, fs); % Coherence: PPG vs. PPGcont

[Cx, Fx] = mscohere(PPGcont, x, window, noverlap, nfft, fs); % Coherence: PPGcont vs. x
[Cy, Fy] = mscohere(PPGcont, y, window, noverlap, nfft, fs); % Coherence: PPGcont vs. y
[Cz, Fz] = mscohere(PPGcont, z, window, noverlap, nfft, fs); % Coherence: PPGcont vs. z
[Cn, Fn] = mscohere(PPGcont, n, window, noverlap, nfft, fs); % Coherence: PPGcont vs. n

figure(4)
plot(F_P_R,C_P_R,'LineWidth',2,'Color',[255/255,144/255,0]); hold on
plot(F_Pcont_R,C_Pcont_R,'LineWidth',2,'Color',[0,206/255,209/255]);
xlabel('Frecuencia (Hz)'); xlim([0 10]); ylim([0 1.1]); grid on
legend('PPG referencia y Recuperada','PPG contaminada y Recuperada'); 
ylabel('Coherencia espectral');

figure(5)
plot(Fx,Cx,'LineWidth',1.7); hold on
plot(Fy,Cy,'LineWidth',1.7); hold on
plot(Fz,Cz,'LineWidth',1.7); hold on
plot(Fn,Cn,'LineWidth',1.7);
xlabel('Frecuencia (Hz)'); xlim([0 10]); grid on
legend('PPG contaminada y x','PPG contaminada y y',...
    'PPG de referencia y z','PPG contaminada y n'); 
ylabel(ax2, 'Coherencia espectral');

figure(6)
plot(PPG,'LineWidth',1.7); hold on
plot(PPGcont+1,'LineWidth',1.7); hold on
plot(R+1,'LineWidth',1.7); grid on
legend({'PPG ref','PPG cont','R'})
xlabel('muestras'); ylabel('A (v)');
%xlim([8300 8900])

%%%%%%%%%%%%%%%% Spectral Error Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency range of interest (sample indices)
a = 26; % lower limit 
b = 37; % upper limit

% Mean power in the selected band
PSDpPPG = sum( PSD_PPGA(a:b) ) / length(PSD_PPGA(a:b));
PSDpPPGcont = sum( PSD_PPGcont(a:b) ) / length(PSD_PPGcont(a:b));
PSDpR = sum( PSD_R(a:b) ) / length(PSD_R(a:b));

% Relative error (%) w.r.t. clean PPG
E1 = ( (PSDpR - PSDpPPG) / PSDpPPG ) * 100;       % Recovered
E2 = ( (PSDpPPGcont - PSDpPPG) / PSDpPPG ) * 100; % Contaminated

%%%%%%%%%%%%%%%%%%%%%% SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snr1 = 10*log10(mean(PPGcont.^2) / mean(ruido.^2));   % Contaminated vs noise
snr2 = 10*log10(mean(RR.^2) / mean(ruido.^2));        % Recovered vs noise

snr3 = 10*log10(mean(PPG.^2) / mean(PPGcont.^2));     % Ref vs contaminated
snr4 = 10*log10(mean(PPG.^2) / mean(RR.^2));          % Ref vs recovered
