
% Carpeta donde se ha descargado fieldtrip
fieldtrip_path='/Users/carolinamigliorelli/Documents/MATLAB/Epilepsy-MEG/libs/fieldtrip';
addpath(fieldtrip_path)
ft_defaults; % Configuración fieldtrip

% Donde están guardados los datos
data_path='E:\HSJD - EEG';
% Nombre del fichero
file_name='EEG_1890.TRC';

hdr = ft_read_header([data_path filesep file_name]);

Label = hdr.label(1:21); % 21 canales estándar
Fs = hdr.Fs;

% Descartamos Oz y Fz
chan_idx=~ismember(Label,'Oz') & ~ismember(Label,'Fpz');
Label(~chan_idx)=[];

% Parámetros para lectura de señal

muestra_crisis=2450; % Muestra del inicio de la crisis
ini = 1.5; % Segundos antes de la muestra
fi = 5; % Segundos después de la muestra

ini_m = muestra_crisis - (1.5 * Fs);
fi_m = muestra_crisis + (5 * Fs);



data = ft_read_data([data_path filesep file_name],'chanindx',chan_idx,'begsample',ini_m,'endsample',fi_m);

overlapping = 0.5; % Solapamiento entre ventanas
ventana = 0.125; % Duración ventana (en segundos) (he puesto 125ms porque es multiplo de la fs y queda mejor)
tiempo = linspace(-ini,fi,size(data,2)); % Vector de tiempo

%% Dibujar señales 
fig=figure;
fig.Color=[1 1 1];
k=max(max(data));
for i = 1:size(data,1)
    plot(tiempo,data(i,:)-k*i,'k')
    hold on
end
axis tight
ax=gca;
ax.YTick=linspace((-i)*k,0,i);
ax.YTickLabel=Label(i:-1:1);
plot([0 0],ax.YLim,'b--')

%% Dibujar mapas topográficos de la suma del valor absoluto

window = round(ventana*Fs);
n=0;
fig=figure;
fig.Color=[1 1 1];
fi=0;

n_plots=floor(((size(data,2)-1)/window - 1)/(1-overlapping));
a=7; b=8; % parámetros subplot
for n=0:n_plots
    ini=(n*(1-overlapping)*window)+1;
    fi=(n*(1-overlapping)+1)*window+1;
    
    topo_val=sum(abs(data(:,ini:fi)),2);
    if n<a*b
        subplot(a,b,n+1)
     elseif n==a*b
        fig=figure;
        fig.Color=[1 1 1];
        subplot(a,b,1)
    else
        subplot(a,b,n-(a*b)+1)
    end
        
    %topoplot(topo_val,'19can.locs','electrodes','labels') % dibujar el mapa con las labels
    topoplot(topo_val,'19can.locs','electrodes','off'); % dibujar el mapa sin las labels ni electrodos
    ax=gca;
    ax.FontSize=8;
    title([num2str(round(tiempo(ini),2)) 's - ' num2str(round(tiempo(fi),2)) 's.'])

end



