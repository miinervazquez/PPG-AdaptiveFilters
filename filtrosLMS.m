PPGo = PPG3; % Select PPG data from one of the available subjects (e.g., PPG1 to PPG7)
PPG = PPGo(:,1); % Reference PPG signal
PPGcont = PPGo(:,2); % Contaminated PPG signal (affected by motion artifacts)
x = PPGo(:,3); % Acceleration signals from X-axis
y = PPGo(:,4); % Acceleration signals from Y-axis
z = PPGo(:,5); % Acceleration signals from Z-axis

% Remove DC offset from all signals
PPG=PPG-mean(PPG);
x=x-mean(x);
y=y-mean(y);
z=z-mean(z);
n = sqrt((x.^2)+(y.^2)+(z.^2)); % Acceleration magnitude from the three axes

%%%%%%%%%%%%%%%%%%% Parametros iniciales del filtro %%%%%%%%%%%%%%%%%%%%%%%
k=length(PPG); 
N=200;
miu = 0.02 ;
gamma = 0.0005;

mov=[x y z n y z]; %Matriz de señales de movimiento: x,y,z,n,xy,xyz
R=zeros(k,6); %Vectores recuperados: [Rx;Ry;Rz;Rn;Rxy;Rxyz]

for j=1:6
    W=zeros(N,1); %Vector de pesos
    e=zeros(1,k); %Vector de error
    P=zeros(N,1);
    r = mov(:,j); %Señal de mov. a utilizar
    
    if (j==5) 
        noisy = R(:,1);
    elseif (j==6) 
        noisy = R(:,5);
    else
        noisy = PPGcont;
    end
    
for i=N:k
    delay = r(i:-1:i-N+1);
    e(i) = noisy(i)-W'*delay; 
    W = lms(e(i),delay,W,miu); %LMS
%     W = slms(e(i),delay,W,miu); %SLMS
%     W = sslms(e(i),delay,W,miu); %Sign-Sign
%     W = nlms(e(i),delay,W,miu,gamma); %NLMS
%     W = vssnlms(e(i),delay,W,miu,C,alpha,P); %VSS NLMS
end
    R(:,j) = e';
end

t = ( 0 : length( PPG ) - 1 ) / 100 ;  %tiempo discreto

figure(1)
plot(t, PPG, 'LineWidth', 1.7); hold on
plot(t, PPGcont + 6, 'LineWidth', 1.7); hold on

% Suponiendo que R tiene las mismas filas que PPG
plot(t, R(:,1) + 1, 'LineWidth', 1.7); hold on
plot(t, R(:,2) + 2, 'LineWidth', 1.7); hold on
plot(t, R(:,3) + 3, 'LineWidth', 1.7); hold on
plot(t, R(:,4) + 4, 'LineWidth', 1.7); hold on
plot(t, R(:,5) + 5, 'LineWidth', 1.7); grid on

legend({'PPG ref', 'PPG cont', 'Rx', 'Ry', 'Rz', 'Rn', 'Rxyz'})
xlabel('Tiempo (s)'); ylabel('Amplitud (v)');
% xlim([82.13 115.92]); %ylim([-0.5 7])

R = R(:,2); %Mejor recuperada
ruido = y ; %Señal de movimiento utilizada

%%%%%%%%%%%%%%%%%%%%%%%%% Espectro de potencias %%%%%%%%%%%%%%%%%%%%%%%%%%%
window=boxcar(128); %Ventana rectangular
noverlap=64; %Solapamiento del 50%
nfft=512; %Tamaño de las secciones

[PSD_PPG, f_PPG]=pwelch(PPG,window,noverlap,nfft,fs);
[PSD_PPGcont, f_PPGcont]=pwelch(PPGcont,window,noverlap,nfft,fs); 
[PSD_R, f_R]=pwelch(RR,window,noverlap,nfft,fs);
[PSD_r, f_r]=pwelch(ruido,window,noverlap,nfft,fs);

PSD_PPGA = PSD_PPG * (max(PSD_PPGcont(4:end))/max(PSD_PPG(4:end))) ;
colores = lines(5); 
espaciado_marcadores = 1;
tamano_marcadores = 10;

figure(2)
ax1 = axes('Position',[0.2 0.15 0.7 0.75]);  
plot(ax1, f_PPG, PSD_PPGA, 'LineWidth', 7, 'Marker', 'o', 'MarkerSize', tamano_marcadores, 'Color', colores(1,:), 'MarkerIndices', 1:espaciado_marcadores:N); hold on
plot(ax1, f_PPGcont, PSD_PPGcont, 'LineWidth', 7, 'Marker', 's', 'MarkerSize', tamano_marcadores, 'Color', colores(2,:), 'MarkerIndices', 1:espaciado_marcadores:N);
plot(ax1, f_R, PSD_R, 'LineWidth', 7, 'Marker', '^', 'MarkerSize', tamano_marcadores, 'Color', colores(3,:), 'MarkerIndices', 1:espaciado_marcadores:N); 
xlim(ax1, [0 10]); 
xlabel(ax1, 'Frecuencia (Hz)', 'FontSize', 40); 
ylabel(ax1, 'PSD (mV^2/Hz)', 'FontSize', 40); 
title(ax1, 'Densidad espectral de potencia', 'FontSize', 50);
grid on
set(ax1, 'Box', 'off', 'TickDir', 'out', 'FontSize', 35);
ax2 = axes('Position',[0.52 0.35 0.35 0.5]);
x_shaded = [3, 5, 5, 3];
y_shaded = [0, 0, 0.03, 0.03];
fill(ax2, x_shaded, y_shaded, [0.9, 0.4, 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % Rojo claro translúcido
hold on;
plot1 = plot(ax2, f_PPG, PSD_PPGA, 'LineWidth', 7, 'Marker', 'o', 'MarkerSize', tamano_marcadores, 'Color', colores(1,:), 'MarkerIndices', 1:espaciado_marcadores:N); 
plot2 = plot(ax2, f_PPGcont, PSD_PPGcont, 'LineWidth', 7, 'Marker', 's', 'MarkerSize', tamano_marcadores, 'Color', colores(2,:), 'MarkerIndices', 1:espaciado_marcadores:N);
plot3 = plot(ax2, f_R, PSD_R, 'LineWidth', 7, 'Marker', '^', 'MarkerSize', tamano_marcadores, 'Color', colores(3,:), 'MarkerIndices', 1:espaciado_marcadores:N); 
grid on
xlim(ax2, [2.34 7]); 
xlabel(ax2, 'Frecuencia (Hz)', 'FontSize', 30); 
ylabel(ax2, 'PSD (mV^2/Hz)', 'FontSize', 30);
legend(ax2, [plot1, plot2, plot3], {'Pulso de referencia', 'Pulso contaminado', 'Pulso recuperado'}, ...
    'FontSize', 30, 'Box', 'off', 'Location', 'north');
set(ax2, 'Box', 'off', 'TickDir', 'out', 'FontSize', 30);

figure(3)
plot6_handle = plot(f_r, PSD_r, 'LineWidth', 7, 'Color', colores(1,:)); 
hold on;
x_shaded_6 = [3, 5, 5, 3]; 
y_shaded_6 = [0, 0, 3e-3, 3e-3]; 
fill(x_shaded_6, y_shaded_6, [0.9, 0.4, 0.7], 'FaceAlpha', 0.1, 'EdgeColor', 'none'); 
xlim([0 10]); 
title('Densidad espectral de potencia', 'FontSize', 10); 
xlabel('Frecuencia (Hz)', 'FontSize', 25); 
ylabel('PSD (mV^2/Hz)', 'FontSize', 30); 
grid on; 
legend('Señal de movimiento en eje y', 'FontSize', 40, 'Box', 'off', 'Location', 'northwest');
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 45);

%%%%%%%%%%%%%%%%%%%%%% Coherencia espectral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PPGcont=PPGcont-mean(PPGcont);
RR = R-mean(R);
fs = 100; %Frecuencia de muestreo

[C_P_R,F_P_R] = mscohere(PPG,RR,window,noverlap,nfft,fs); %Coherencia PPG y Recuperada
[C_Pcont_R,F_Pcont_R] = mscohere(PPGcont,RR,window,noverlap,nfft,fs); %Coherencia PPGcont y Recuperada
[C_P_Pcont,F_P_Pcont] = mscohere(PPG,PPGcont,window,noverlap,nfft,fs); %Coherencia PPG y PPGcont

[Cx,Fx] = mscohere(PPGcont,x,window,noverlap,nfft,fs); %Coherencia PPGcont y x
[Cy,Fy] = mscohere(PPGcont,y,window,noverlap,nfft,fs); %Coherencia PPGcont y y
[Cz,Fz] = mscohere(PPGcont,z,window,noverlap,nfft,fs); %Coherencia PPGcont y z
[Cn,Fn] = mscohere(PPGcont,n,window,noverlap,nfft,fs); %Coherencia PPGcont y n

figure(4)
plot(F_P_R,C_P_R,'LineWidth',2,'Color',[255/255,144/255,0]); hold on
plot(F_Pcont_R,C_Pcont_R,'LineWidth',2,'Color',[0,206/255,209/255]);
xlabel('Frecuencia (Hz)'); xlim([0 10]); ylim([0 1.1]); grid on
legend('PPG referencia y Recuperada','PPG contaminada y Recuperada'); 
ylabel('Coherencia espectral');

figure(5)
ax1 = axes('Position',[0.2 0.1 0.7 0.8]);
plot(ax1, Fx,Cx,Fy,Cy,Fz,Cz,Fn,Cn,'LineWidth',2); 
xlim(ax1, [0 10]); % ylim(ax1, [0 0.0045]);
legend(ax1, 'PPG contaminada y x','PPG contaminada y y',...
    'PPG de referencia y z','PPG contaminada y n'); grid on
xlabel(ax1, 'Frecuencia (Hz)'); ylabel(ax1, 'Coherencia espectral');
%%%%
ax2 = axes('Position',[0.5 0.35 0.38 0.4]); %[0.64 0.28 0.25 0.5]
plot(ax2, Fx,Cx,'LineWidth',2); hold on; 
plot(ax2, Fy,Cy,'LineWidth',2); hold on
plot(ax2, Fz,Cz,'LineWidth',2); grid on
plot(ax2, Fn,Cn,'LineWidth',2); grid on
xlim(ax2, [5 7]); % ylim(ax2, [0 0.006]);
xlabel(ax2, 'Frecuencia (Hz)'); ylabel(ax2, 'Coherencia espectral'); 

figure(6)
plot(PPG,'LineWidth',1.7); hold on
plot(PPGcont+1,'LineWidth',1.7); hold on
plot(R+1,'LineWidth',1.7); grid on
legend({'PPG ref','PPG cont','R'})
xlabel('muestras'); ylabel('A (v)');
xlim([8300 8900])

%%%%%%%%%%%%%%%% Calculo del error espectral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rango de frecuencias filtradas - no. de muestra correspondiente
a = 26; %limite inferior 
b = 37; %limite superior

PSDpPPG = sum( PSD_PPGA(a:b) ) / length(PSD_PPGA(a:b));
PSDpPPGcont = sum( PSD_PPGcont(a:b) ) / length(PSD_PPGcont(a:b));
PSDpR = sum( PSD_R(a:b) ) / length(PSD_R(a:b));

E1 = ( (PSDpR - PSDpPPG) / PSDpPPG ) * 100 ;
E2 = ( (PSDpPPGcont - PSDpPPG) / PSDpPPG ) * 100 ;

%%%%%%%%%%%%%%%%%%%%%% SNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snr1=10*log10(mean(PPGcont.^2)/mean(ruido.^2));
snr2=10*log10(mean(RR.^2)/mean(ruido.^2));

snr3=10*log10(mean(PPG.^2)/mean(PPGcont.^2));
snr4=10*log10(mean(PPG.^2)/mean(RR.^2));

