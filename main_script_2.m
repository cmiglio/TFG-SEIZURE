close all
clear all
clc

% Carpeta donde se ha descargado fieldtrip
fieldtrip_path='/Users/beatrizalmajano/Documents/fieldtrip-master';
addpath(fieldtrip_path)
ft_defaults; % Configuración fieldtrip

% Donde están guardados los datos
data_path='/Volumes/Elements SE/TFG-epilepsia/PAT_7';
% Nombre del fichero
file_name='EEG_1050.TRC';

hdr = ft_read_header([data_path filesep file_name]);

Label = hdr.label(1:21); % 21 canales estándar
Fs = hdr.Fs;

% Descartamos Oz y Fz
chan_idx=~ismember(Label,'Oz') & ~ismember(Label,'Fpz');
Label(~chan_idx)=[];

% Parámetros para lectura de señal
segundos_crisis=5908.652;
muestra_crisis=round(segundos_crisis*Fs); % Muestra del inicio de la crisis
ini = 1.5; % Segundos antes de la muestra
fi = 5; % Segundos después de la muestra

ini_m = muestra_crisis - (1.5 * Fs);
fi_m = muestra_crisis + (5 * Fs);


win_filt=128;
data = ft_read_data([data_path filesep file_name],'chanindx',chan_idx,'begsample',ini_m-win_filt,'endsample',fi_m+win_filt);

Hd = filtro_lp_40(Fs);

data=filtfilt(Hd.sosMatrix,Hd.ScaleValues,data');
data=data(win_filt+1:end-128,:)';


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

overlapping=0;
window=15; % 31, ms
n_plots=floor(((size(data,2)-1)/window - 1)/(1-overlapping));


fig=figure;
for n=1:n_plots
    ini=(n*(1-overlapping)*window)+1;
    fi=(n*(1-overlapping)+1)*window+1;
    
    %% Dibujar señales 
%         subplot(1,2,1)
%         k=max(max(data));
%         for i = 1:size(data,1)
%             plot(tiempo,data(i,:)-k*i,'k')
%             hold on
%         end
%         axis tight
%         ax=gca;
%         ax.YTick=linspace((-i)*k,0,i);
%         ax.YTickLabel=Label(i:-1:1);
%         plot([tiempo(ini) tiempo(ini)],ax.YLim,'b--')
%         hold off
    
    
    topo_val=sum(abs(data(:,ini:fi)),2);
    topo_val2(:,n)=sum(abs(data(:,ini:fi)),2);
        
    %topoplot(topo_val,'19can.locs','electrodes','labels') % dibujar el mapa con las labels
%     subplot(1,2,2)
%     topoplot(topo_val,'19can.locs','electrodes','on','numcontour',4); % dibujar el mapa sin las labels ni electrodos
%     caxis([0 3000])
%     colorbar('South')
%     ax=gca;
%     ax.FontSize=8;
%     title(num2str(round(tiempo(ini),2)))
    %waitforbuttonpress
end
% Canales frontales
figure;
subplot(3,1,1); plot(topo_val2(1:7,:)')
legend(Label(1:7))
ylim([0 3000])

% Canales centrales
%figure;
subplot(3,1,2); plot(topo_val2(8:12,:)')
legend(Label(8:12))
ylim([0 3000])

% Canales parietales/occipitales
%figure;
subplot(3,1,3); plot(topo_val2(13:19,:)')
legend(Label(13:19))
ylim([0 3000])



