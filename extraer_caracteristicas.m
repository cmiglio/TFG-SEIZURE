close all
clear all
warning off
clc

% Carpeta donde se ha descargado fieldtrip
% fieldtrip_path='/Users/carolinamigliorelli/Documents/MATLAB/Epilepsy-MEG/libs/fieldtrip';
% fieldtrip_path='C:\Users\cmigliorelli\Documents\MATLAB\fieldtrip';
fieldtrip_path='C:\Users\cmigliorelli\Documents\MATLAB\fieldtrip';
addpath(fieldtrip_path)
addpath('Functions')

ft_defaults; % Configuración fieldtrip

% Donde están guardados los datos
%data_path='/Volumes/KINGSTON/Epilepsia_HSJD/PAT_3/EEG';
%data_path='Y:\Electrofisiologia\PAT_3';
data_path='Y:\Electrofisiologia\PAT_7\';
% Nombre del fichero
file_name='EEG_1050.TRC';

% Cargar fichero de eventos
load('C:\Users\cmigliorelli\Documents\MATLAB\seizureDetection\eventos\PAT_7\events_EEG_1050.mat');


hdr = ft_read_header([data_path filesep file_name]);

Label = hdr.label(1:21); % 21 canales estándar
Fs = hdr.Fs;

% Descartamos Oz y Fz
chan_idx=~ismember(Label,'Oz') & ~ismember(Label,'Fpz');
Label(~chan_idx)=[];


cfg = [];
cfg.dataset     = [data_path filesep file_name];
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 0.5;
cfg.dftfilter   = 'yes';
cfg.channel     = Label;
data_org        = ft_preprocessing(cfg);

data            = data_org.trial{1};
time            = data_org.time{1};
es_crisis       = zeros(size(time));


for i=1:size(events.samples,2)
    es_crisis(events.samples(1,i):events.samples(2,i))=1;
end


data = data(:,idx:end);
time = time(idx:end);
es_crisis = es_crisis(idx:end);

% 
ov = 0.5; % Solapamiento entre ventanas
ventana = 2; % Duración ventana (en segundos) (he puesto 125ms porque es multiplo de la fs y queda mejor)
ventana_m = ventana*Fs;
l_total=ceil(size(data,2)/(ventana_m*(1-ov)))-1;


bandas=[0.5 30; 0.5 4; 4 8; 8 12; 12 25];

features=zeros(l_total,size(data,1),3,size(bandas,1)); % 7 parámetros, 5 bandas, 19 canales
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

            [Pxx,frecs] = pwelch(data(j,ini:fin),2^(nextpow2(fin-ini)-1),0.9,2^nextpow2(fin-ini),Fs);

            for k = 1:size(bandas,1) % 5 bandas
                frec_band=bandas(k,:);
                params = calc_params2(Pxx,frecs,Fs,frec_band);
                for n=1:3
                    features(i,j,n,k)=params(n);
                end
            end

        end
        if rem(i,100)==0
            disp(['Completado ' num2str(100*i/l_total) '%'])
        end
    end
end
  
es_crisis_win=int8(es_crisis_win);


% Explore features:
for i = 1:l_total
   topoplot(features(i,:,3,1),'19can.locs','noplot','off');
   colorbar
    title(['crisis: ' num2str(es_crisis_win(i))])
    waitforbuttonpress
end

