close all
clear all
warning off
clc

% Carpeta donde se ha descargado fieldtrip
% fieldtrip_path='/Users/carolinamigliorelli/Documents/MATLAB/Epilepsy-MEG/libs/fieldtrip';
 fieldtrip_path='C:\Users\cmigliorelli\Documents\MATLAB\fieldtrip';
%fieldtrip_path='C:\Users\Caro\Documents\MATLAB\Epilepsy-MEG\libs\fieldtrip';
addpath(fieldtrip_path)
addpath('Functions')

ft_defaults; % Configuración fieldtrip

% Donde están guardados los datos
%data_path='/Volumes/KINGSTON/Epilepsia_HSJD/PAT_3/EEG';
data_path='Y:\Electrofisiologia\PAT_3';
%data_path='E:\Epilepsia_HSJD\PAT_3\EEG';
% Nombre del fichero
file_name='EEG_365.TRC';

results_dir='PAT_3';
mkdir(results_dir)

% Cargar fichero de eventos
load('events_MarkerFile-bst.mat');


hdr = ft_read_header([data_path filesep file_name]);

Label = hdr.label(1:21); % 21 canales estándar
Fs = hdr.Fs;

% Descartamos Oz y Fz
chan_idx=~ismember(Label,'Oz') & ~ismember(Label,'Fpz');
Label(~chan_idx)=[];


cfg = [];
cfg.dataset     = [data_path filesep file_name];
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [0.5 30];
cfg.channel     = Label;
data_org        = ft_preprocessing(cfg);

data            = data_org.trial{1};
time            = data_org.time{1};
es_crisis       = zeros(size(time));


for i=1:size(events.samples,2)
    es_crisis(events.samples(1,i):events.samples(2,i))=1;
end

idx=find(es_crisis==1,1,'first');

data = data(:,idx:end);
time = time(idx:end);
es_crisis = es_crisis(idx:end);

% 
ov = 0; % Solapamiento entre ventanas
ventana = 2; % Duración ventana (en segundos) (he puesto 125ms porque es multiplo de la fs y queda mejor)
ventana_m = ventana*Fs;
l_total=ceil(size(data,2)/(ventana_m*(1-ov)))-1;


bandas=[0.5 30; 0.5 4; 4 8; 8 12; 12 25];

features=zeros(l_total,size(data,1),8); % 8 parámetros, 19 canales
time_win =zeros(1,l_total);
 es_crisis_win=zeros(1,l_total);

% Features:

% 1. Median frequency
% 2. RMS
% 3. Sample entropy


for i = 1:l_total
    ini=(i-1)*(1-ov)*ventana_m+1;
    fin=ini+ventana_m;
    
    if fin<=length(time)
        time_win(i) = mean(time(ini:fin));  
        es_crisis_win(i) =  mean(es_crisis(ini:fin))>=0.5;
        
        for j=1:size(data,1)   

            X=data(j,ini:fin);
            [Pxx,f] = periodogram(X,[],[],Fs);
            for k = 1:size(bandas,1) % 5 bandas
                features(i,j,k) = bandpower(X, Fs, bandas(k,:));
            end
                features(i,j,k+1)=SampEn(2 ,0.2*std(X),X);
                features(i,j,k+2)=rms(X);
                features(i,j,k+3)=medfreq(Pxx,f);
        end
        if rem(i,100)==0
            disp(['Completado ' num2str(100*i/l_total) '%'])
        end
    end
end
  
es_crisis_win=int8(es_crisis_win);


% Valores medios durante la crisis
power_cr=mean(features(es_crisis_win==1,:,1));
power_delta_cr=mean(features(es_crisis_win==1,:,2));
power_theta_cr=mean(features(es_crisis_win==1,:,3));
power_alpha_cr=mean(features(es_crisis_win==1,:,4));
power_beta_cr=mean(features(es_crisis_win==1,:,5));
entropy_cr=mean(features(es_crisis_win==1,:,6));
rms_cr=mean(features(es_crisis_win==1,:,7));
mfrec_cr=mean(features(es_crisis_win==1,:,8));

% Baseline
power_tot=mean(features(:,:,1));
power_delta_tot=mean(features(:,:,2));
power_theta_tot=mean(features(:,:,3));
power_alpha_tot=mean(features(:,:,4));
power_beta_tot=mean(features(:,:,5));
entropy_tot=mean(features(:,:,6));
rms_tot=mean(features(:,:,7));
mfrec_tot=mean(features(:,:,8));

% No crisis
power_ncr=mean(features(es_crisis_win==0,:,1));
power_delta_ncr=mean(features(es_crisis_win==0,:,2));
power_theta_ncr=mean(features(es_crisis_win==0,:,3));
power_alpha_ncr=mean(features(es_crisis_win==0,:,4));
power_beta_ncr=mean(features(es_crisis_win==0,:,5));
entropy_ncr=mean(features(es_crisis_win==0,:,6));
rms_ncr=mean(features(es_crisis_win==0,:,7));
mfrec_ncr=mean(features(es_crisis_win==0,:,8));


save([results_dir filesep 'features.mat'],'features')

% Topografias durante crisis
fig1=figure(1);
subplot(2,4,1)
topoplot(power_cr./power_tot,'19can.locs');
title('Potencia')
colorbar
caxis([0 max(power_cr./power_tot)])
subplot(2,4,2)
topoplot(power_delta_cr./power_delta_tot,'19can.locs');
title('Potencia delta total')
colorbar
caxis([0 max(power_delta_cr./power_delta_tot)])
subplot(2,4,3)
topoplot(power_theta_cr./power_theta_tot,'19can.locs');
title('Potencia theta total')
colorbar
caxis([0 max(power_theta_cr./power_theta_tot)])
subplot(2,4,4)
topoplot(power_alpha_cr./power_alpha_tot,'19can.locs');
caxis([0 max(power_alpha_cr./power_alpha_tot)])
title('Potencia alpha total')
colorbar
subplot(2,4,5)
topoplot(power_beta_cr./power_beta_tot,'19can.locs');
title('Potencia beta total')
colorbar
caxis([0 max(power_beta_cr./power_beta_tot)])
subplot(2,4,6)
topoplot(entropy_cr./entropy_tot,'19can.locs');
title('Sample entropy')
colorbar
caxis([min(entropy_cr./entropy_tot) max(entropy_cr./entropy_tot)])
subplot(2,4,7)
topoplot(rms_cr./rms_tot,'19can.locs');
title('RMS')
colorbar
caxis([0 max(rms_cr./rms_tot)])
subplot(2,4,8)
topoplot(mfrec_cr./mfrec_tot,'19can.locs');
title('Frecuencia mediana')
colorbar
caxis([0 max(mfrec_cr./mfrec_tot)])

saveas(fig1,[results_dir filesep 'fig_crisis.png'])

% Topografias no crisis
fig2=figure(2);
subplot(2,4,1)
topoplot(power_ncr./power_tot,'19can.locs');
title('Potencia')
colorbar
caxis([0 max(power_ncr./power_tot)])
subplot(2,4,2)
topoplot(power_delta_ncr./power_delta_tot,'19can.locs');
title('Potencia delta total')
colorbar
caxis([0 max(power_delta_ncr./power_delta_tot)])
subplot(2,4,3)
topoplot(power_theta_ncr./power_theta_tot,'19can.locs');
title('Potencia theta total')
colorbar
caxis([0 max(power_theta_ncr./power_theta_tot)])
subplot(2,4,4)
topoplot(power_alpha_ncr./power_alpha_tot,'19can.locs');
caxis([0 max(power_alpha_ncr./power_alpha_tot)])
title('Potencia alpha total')
colorbar
subplot(2,4,5)
topoplot(power_beta_ncr./power_beta_tot,'19can.locs');
title('Potencia beta total')
colorbar
caxis([0 max(power_beta_ncr./power_beta_tot)])
subplot(2,4,6)
topoplot(entropy_ncr./entropy_tot,'19can.locs');
title('Sample entropy')
colorbar
caxis([min(entropy_ncr./entropy_tot) max(entropy_ncr./entropy_tot)])
subplot(2,4,7)
topoplot(rms_ncr./rms_tot,'19can.locs');
title('RMS')
colorbar
caxis([0 max(rms_ncr./rms_tot)])
subplot(2,4,8)
topoplot(mfrec_ncr./mfrec_tot,'19can.locs');
title('Frecuencia mediana')
colorbar
caxis([0 max(mfrec_ncr./mfrec_tot)])


saveas(fig2,[results_dir filesep 'fig_nocrisis.png'])

