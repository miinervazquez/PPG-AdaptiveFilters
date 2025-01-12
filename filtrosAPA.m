
PPG=PPG-mean(PPG);
PPGcont=PPGcont-mean(PPGcont);
x=x-mean(x);
y=y-mean(y);
z=z-mean(z);
n = sqrt((x.^2)+(y.^2)+(z.^2));

%%%%%%%%%%%%%%%%%%% Parametros iniciales del filtro %%%%%%%%%%%%%%%%%%%%%%%
k = length(PPG); 
N = 200;
M = 10;

miu = 0.03; % mu2max is chosen less than 2
gamma = 0.002;
C = 0.002 ; % C is a positive constant.
alpha = 0.2 ; % 0 < alpha < 1

mov=[x y z n y z]; %Matriz de señales de movimiento: x,y,z,n,xy,xyz
R=zeros(k,6); %Vectores recuperados: [Rx;Ry;Rz;Rn;Rxy;Rxyz]

for j=1:6
    W=zeros(N,1); %Vector de pesos
    e=zeros(1,k); %Vector de error
    X=zeros(N,M);
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
    delay=r(i:-1:i-N+1);
    X = [delay X(:,1:M-1)];
    Y = X'*W;
    E = noisy(i:-1:i-M+1)-Y; 
    e(i)=mean(E);
    W = apa(W,miu,X,E,M,gamma); %APA
%     W = vapa(W,miu,X,E); %VAPA
%     W = vssapa(W,miu,X,E,C,alpha,P,M,gamma); %Variable Step Size APA
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
plot(t, R(:,6) + 5, 'LineWidth', 1.7); grid on

legend({'PPG ref', 'PPG cont', 'Rx', 'Ry', 'Rz', 'Rn', 'Rxyz'})
xlabel('Tiempo (s)'); ylabel('Amplitud (v)');
% xlim([82.13 115.92]); %ylim([-0.5 7])

R = R(:,2); %Mejor recuperada
ruido = y ; %Señal de movimiento utilizada
RR = R-mean(R);

%%%%%%%%%%%%%%%%%%%%%%%%% Espectro de potencias %%%%%%%%%%%%%%%%%%%%%%%%%%% 
window=boxcar(128); %Ventana rectangular
noverlap=64; %Solapamiento del 50%
nfft=512; %Tamaño de las secciones
fs = 100; %Frecuencia de muestreo

[PSD_PPG, f_PPG]=pwelch(PPG,window,noverlap,nfft,fs);
[PSD_PPGcont, f_PPGcont]=pwelch(PPGcont,window,noverlap,nfft,fs); 
[PSD_R, f_R]=pwelch(RR,window,noverlap,nfft,fs);
[PSD_r, f_r]=pwelch(ruido,window,noverlap,nfft,fs);

PSD_PPGA = PSD_PPG * (max(PSD_PPGcont(4:end))/max(PSD_PPG(4:end))) ;

figure(2)
ax1 = axes('Position',[0.2 0.1 0.7 0.8]);
plot(ax1, f_PPG,PSD_PPGA,f_PPGcont,PSD_PPGcont,f_R,PSD_R,'LineWidth',2); 
xlim(ax1, [0 10]); % ylim(ax1, [0 0.0045]);
legend(ax1, 'Pulso de referencia', 'Pulso contaminado','Pulso recuperado')
xlabel(ax1, 'Frecuencia (Hz)'); ylabel(ax1, 'PSD (mV^2/Hz)'); grid on
title(ax1, 'Densidad espectral de potencia')
%%%%
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

%%%%%%%%%%%%%%%%%%%%%% Coherencia espectral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%% Calculo del error espectral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rango de frecuencias filtradas - no. de muestra correspondiente
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

