close all
clear all
warning off
clc


% Carpeta donde se ha descargado fieldtrip
% fieldtrip_path='/Users/carolinamigliorelli/Documents/MATLAB/Epilepsy-MEG/libs/fieldtrip';
%fieldtrip_path='C:\Users\cmigliorelli\Documents\MATLAB\fieldtrip';
%fieldtrip_path='D:\Programas\fieldtrip-20180213';
fieldtrip_path='C:\Users\Caro\Documents\MATLAB\Epilepsy-MEG\libs\fieldtrip';
addpath(fieldtrip_path)
addpath('Functions')

ft_defaults; % Configuración fieldtrip


patient={'PAT_7','PAT_7','PAT_9','PAT_9','PAT_12','PAT_12','PAT_13','PAT_13'};
filenames={'EEG_1059','EEG_1050','EEG_1275','EEG_1220','EEG_1692','EEG_1681','EEG_1933','EEG_1923'};

for pat=1:length(patient)

    % Donde están guardados los datos
    %data_path='/Volumes/KINGSTON/Epilepsia_HSJD/PAT_3/EEG';
    data_path=['Y:\Electrofisiologia' filesep patient{pat}];
    data_path=['D:\BBDD\BIOART_signals\Epilepsia_HSJD\BBDD\trc_data' filesep patient{pat}];
    %data_path='E:\Epilepsia_HSJD\PAT_3\EEG';
    % Nombre del fichero
    file_name=filenames{pat};
    
    results_dir=['results' patient{pat} '-' filenames{pat}];
    mkdir(results_dir)

    % Cargar fichero de eventos
    load(['eventos' filesep patient{pat} filesep 'events_' filenames{pat} '.mat']);


    hdr = ft_read_header([data_path filesep file_name '.TRC']);

    Label = hdr.label(1:21); % 21 canales estándar
    Fs = hdr.Fs;

    % Descartamos Oz y Fz
    chan_idx=~ismember(Label,'Oz') & ~ismember(Label,'Fpz');
    Label(~chan_idx)=[];


    cfg = [];
    cfg.dataset     = [data_path filesep file_name '.TRC'];
    %cfg.lpfilter    = 'yes';
    %cfg.lpfreq      = 0.5;
    cfg.dftfreq     = [50, 100, 150];
    cfg.dftfilter    = 'yes';
    cfg.channel     = Label;
    data_org        = ft_preprocessing(cfg);

    data            = data_org.trial{1};
    time            = data_org.time{1};
    es_crisis       = zeros(size(time));


    % Se cargan las features ya calculadas. SOLAMENTE QUIERO CAMBIAR LA
    % VARIABLE "es_crisis_win"
    load(['results' patient{pat} '-' filenames{pat} filesep 'features.mat'])
    
    
    
    for i=1:size(events.samples,2)
        es_crisis(events.samples(1,i):events.samples(2,i))=1;        
    end


    
% %     for i=1:size(events.samples,2)
% %         
% %         ini=max(1,events.samples(1,i)-40*Fs);
% %         fin=min(size(es_crisis,2),events.samples(2,i)+40*Fs);
% %         if sum(es_crisis(ini:events.samples(1,i)-1))==0
% %             es_crisis(ini:events.samples(1,i)-1)=1; 
% %         end 
% %         if sum(es_crisis(events.samples(2,i)+1:fin))==0
% %             es_crisis(events.samples(2,i)+1:fin)=3; 
% %         end
% %         
% %     end

 
    ov = 0; % Solapamiento entre ventanas
    ventana = 5; % Duración ventana (en segundos) 
    ventana_m = ventana*Fs;
    l_total=ceil(size(data,2)/(ventana_m*(1-ov)))-1;

    % Añadir las bandas de antes también. 
    bandas=[0.5 30; 0.5 4; 4 8; 8 12; 12 25; 25 45; 55 80;80 95; 105 120];
    % gamma baja(25-45), Gamma alta(55 80), ripple1(80 95), ripple2(105 128)
    %features=zeros(l_total,size(data,1),12); % 12 parámetros, 19 canales
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
            es_crisis_win(i) =  round(mean(es_crisis(ini:fin)));         
            
% %             for j=1:size(data,1)   
% % 
% %                 X=data(j,ini:fin);
% %                 [Pxx,f] = periodogram(X,[],[],Fs);
% %                 for k = 1:size(bandas,1) % 9 bandas 
% %                     features(i,j,k) = mean(Pxx(f>bandas(k,1) & f<bandas(k,2)));
% %                 end
% %                     features(i,j,k+1)=SampEn(2 ,0.2*std(X),X);
% %                     features(i,j,k+2)=rms(X);
% %                     features(i,j,k+3)=medfreq(Pxx,f);
% %             end
            if rem(i,100)==0
                disp(['Completado ' num2str(100*i/l_total) '%'])
            end
        end
    end
    
    
    
    [~,win_ini] = find(diff(es_crisis_win)==1);
    [~,win_fin] = find(diff(es_crisis_win)==-1);
    
    for n=1:size(win_ini,2)
        if (sum(es_crisis_win(win_ini(n)-7:win_ini(n)))==0)
            es_crisis_win(win_ini(n)-7:win_ini(n))=2; 
        end 

        if (sum(es_crisis_win(win_fin(n)+1:win_fin(n)+8))==0)
            es_crisis_win(win_fin(n)+1:win_fin(n)+8)=3; 
        end 
    end
    
    
      
   
    save([results_dir filesep 'features2.mat'],'features','es_crisis_win')

  
end
