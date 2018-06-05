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


patient={'PAT_7','PAT_7','PAT_9','PAT_9','PAT_12','PAT_12','PAT_13','PAT_13'};
filenames={'EEG_1059','EEG_1050','EEG_1275','EEG_1220','EEG_1692','EEG_1681','EEG_1933','EEG_1923'};

for pat=1:length(patient)

    % Donde están guardados los datos
    %data_path='/Volumes/KINGSTON/Epilepsia_HSJD/PAT_3/EEG';
    data_path=['Y:\Electrofisiologia' filesep patient{pat}];
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
    %cfg.dftfilter    = 'yes';
    cfg.channel     = Label;
    data_org        = ft_preprocessing(cfg);

    data            = data_org.trial{1};
    time            = data_org.time{1};
    es_crisis       = zeros(size(time));


    for i=1:size(events.samples,2)
        
        es_crisis(events.samples(1,i):events.samples(2,i))=1;
        ini=max(1,events.samples(1,i)-40*Fs);
        fin=min(size(es_crisis,2),events.samples(2,i)+40*Fs);
        if sum(es_crisis(ini:events.samples(1,i)-1))==0
            es_crisis(ini:events.samples(1,i)-1)=2; 
        end 
        if sum(es_crisis(events.samples(2,i)+1:fin))==0
            es_crisis(events.samples(2,i)+1:fin)=3; 
        end
        
    end

    % 
    ov = 0; % Solapamiento entre ventanas
    ventana = 5; % Duración ventana (en segundos) 
    ventana_m = ventana*Fs;
    l_total=ceil(size(data,2)/(ventana_m*(1-ov)))-1;

    % Añadir las bandas de antes también. 
    bandas=[0.5 30; 0.5 4; 4 8; 8 12; 12 25; 25 45; 55 80;80 95; 105 120];
    % gamma baja(25-45), Gamma alta(55 80), ripple1(80 95), ripple2(105 128)
    features=zeros(l_total,size(data,1),12); % 12 parámetros, 19 canales
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

            for j=1:size(data,1)   

                X=data(j,ini:fin);
                [Pxx,f] = periodogram(X,[],[],Fs);
                for k = 1:size(bandas,1) % 9 bandas 
                    features(i,j,k) = mean(Pxx(f>bandas(k,1) & f<bandas(k,2)));
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

    
    % % Valores medios durante la crisis
    % power_cr=mean(features(es_crisis_win==1,:,1));
    % power_delta_cr=mean(features(es_crisis_win==1,:,2));
    % power_theta_cr=mean(features(es_crisis_win==1,:,3));
    % power_alpha_cr=mean(features(es_crisis_win==1,:,4));
    % power_beta_cr=mean(features(es_crisis_win==1,:,5));
    % power_gammabaja_cr=mean(features(es_crisis_win==1,:,6));
    % power_gammaalta_cr=mean(features(es_crisis_win==1,:,7));
    % power_ripple1_cr=mean(features(es_crisis_win==1,:,8));
    % power_ripple2_cr=mean(features(es_crisis_win==1,:,9));
    % entropy_cr=mean(features(es_crisis_win==1,:,10));
    % rms_cr=mean(features(es_crisis_win==1,:,11));
    % mfrec_cr=mean(features(es_crisis_win==1,:,12));
    % 
    % % Baseline
    % power_tot=mean(features(:,:,1));
    % power_delta_tot=mean(features(:,:,2));
    % power_theta_tot=mean(features(:,:,3));
    % power_alpha_tot=mean(features(:,:,4));
    % power_beta_tot=mean(features(:,:,5));
    % power_gammabaja_tot=mean(features(:,:,6));
    % power_gammaalta_tot=mean(features(:,:,7));
    % power_ripple1_tot=mean(features(:,:,8));
    % power_ripple2_tot=mean(features(:,:,9));
    % entropy_tot=mean(features(:,:,10));
    % rms_tot=mean(features(:,:,11));
    % mfrec_tot=mean(features(:,:,12));
    % 
    % % No crisis
    % power_ncr=mean(features(es_crisis_win==0,:,1));
    % power_delta_ncr=mean(features(es_crisis_win==0,:,2));
    % power_theta_ncr=mean(features(es_crisis_win==0,:,3));
    % power_alpha_ncr=mean(features(es_crisis_win==0,:,4));
    % power_beta_ncr=mean(features(es_crisis_win==0,:,5));
    % power_gammabaja_ncr=mean(features(es_crisis_win==0,:,6));
    % power_gammaalta_ncr=mean(features(es_crisis_win==0,:,7));
    % power_ripple1_ncr=mean(features(es_crisis_win==0,:,8));
    % power_ripple2_ncr=mean(features(es_crisis_win==0,:,9));
    % entropy_ncr=mean(features(es_crisis_win==0,:,10));
    % rms_ncr=mean(features(es_crisis_win==0,:,11));
    % mfrec_ncr=mean(features(es_crisis_win==0,:,12));
    % 

    save([results_dir filesep 'features.mat'],'features','es_crisis_win')

    % Topografias durante crisis (potencia total, RMS, median frequency,sample entropy)

    % fig1=figure(1);
    % 
    % subplot(2,2,1)
    % topoplot(power_cr./power_tot,'19can.locs');
    % title('Potencia')
    % colorbar
    % caxis([0 max(power_cr./power_tot)])
    % 
    % subplot(2,2,2)
    % topoplot(entropy_cr./entropy_tot,'19can.locs');
    % title('Sample entropy')
    % colorbar
    % caxis([min(entropy_cr./entropy_tot) max(entropy_cr./entropy_tot)])
    % 
    % subplot(2,2,3)
    % topoplot(rms_cr./rms_tot,'19can.locs');
    % title('RMS')
    % colorbar
    % caxis([0 max(rms_cr./rms_tot)])
    % 
    % subplot(2,2,4)
    % topoplot(mfrec_cr./mfrec_tot,'19can.locs');
    % title('Frecuencia mediana')
    % colorbar
    % caxis([0 max(mfrec_cr./mfrec_tot)])
    % 
    % saveas(fig1,[results_dir filesep 'fig_crisis1.png'])
    % 
    % % Topografias durante crisis (delta,theta,alpha,beta, gamma baja, gamma alta, ripple 1, ripple 2)
    % 
    % fig2=figure(2);
    % 
    % subplot(2,4,1)
    % topoplot(power_delta_cr./power_delta_tot,'19can.locs');
    % title('Potencia delta total')
    % colorbar
    % caxis([0 max(power_delta_cr./power_delta_tot)])
    % 
    % subplot(2,4,2)
    % topoplot(power_theta_cr./power_theta_tot,'19can.locs');
    % title('Potencia theta total')
    % colorbar
    % caxis([0 max(power_theta_cr./power_theta_tot)])
    % 
    % subplot(2,4,3)
    % topoplot(power_alpha_cr./power_alpha_tot,'19can.locs');
    % title('Potencia alpha total')
    % colorbar
    % caxis([0 max(power_alpha_cr./power_alpha_tot)])
    % 
    % subplot(2,4,4)
    % topoplot(power_beta_cr./power_beta_tot,'19can.locs');
    % title('Potencia beta total')
    % colorbar
    % caxis([0 max(power_beta_cr./power_beta_tot)])
    % 
    % subplot(2,4,5)
    % topoplot(power_gammabaja_cr./power_gammabaja_tot,'19can.locs');
    % title('Potencia Gamma Baja')
    % colorbar
    % caxis([0 max(power_gammabaja_cr./power_gammabaja_tot)])
    % 
    % subplot(2,4,6)
    % topoplot(power_gammaalta_cr./power_gammaalta_tot,'19can.locs');
    % title('Potencia Gamma Alta')
    % colorbar
    % caxis([0 max(power_gammaalta_cr./power_gammaalta_tot)])
    % 
    % subplot(2,4,7)
    % topoplot(power_ripple1_cr./power_ripple1_tot,'19can.locs');
    % title('Potencia Ripple1')
    % colorbar
    % caxis([0 max(power_ripple1_cr./power_ripple1_tot)])
    % 
    % subplot(2,4,8)
    % topoplot(power_ripple2_cr./power_ripple2_tot,'19can.locs');
    % title('Potencia Ripple 2')
    % colorbar
    % caxis([0 max(power_ripple2_cr./power_ripple2_tot)])
    % 
    % saveas(fig2,[results_dir filesep 'fig_crisis2.png'])
    % 
    % % Topografias no crisis (potencia total, RMS, median frequency,sample entropy)
    % 
    % fig3=figure(3);
    % 
    % subplot(2,2,1)
    % topoplot(power_ncr./power_tot,'19can.locs');
    % title('Potencia')
    % colorbar
    % 
    % subplot(2,2,2)
    % topoplot(entropy_ncr./entropy_tot,'19can.locs');
    % title('Sample entropy')
    % colorbar
    % caxis([min(entropy_ncr./entropy_tot) max(entropy_ncr./entropy_tot)])
    % 
    % subplot(2,2,3)
    % topoplot(rms_ncr./rms_tot,'19can.locs');
    % title('RMS')
    % colorbar
    % caxis([0 max(rms_ncr./rms_tot)])
    % 
    % subplot(2,2,4)
    % topoplot(mfrec_ncr./mfrec_tot,'19can.locs');
    % title('Frecuencia mediana')
    % colorbar
    % caxis([0 max(mfrec_cr./mfrec_tot)])
    % 
    % 
    % saveas(fig3,[results_dir filesep 'fig_nocrisis1.png'])
    % 
    % % Topografias no crisis (delta,theta,alpha,beta, gamma baja, gamma alta, ripple 1, ripple 2)
    % 
    % fig4=figure(4);
    % 
    % subplot(2,4,1)
    % topoplot(power_delta_ncr./power_delta_tot,'19can.locs');
    % title('Potencia delta total')
    % colorbar
    % caxis([0 max(power_delta_ncr./power_delta_tot)])
    % 
    % subplot(2,4,2)
    % topoplot(power_theta_ncr./power_theta_tot,'19can.locs');
    % title('Potencia theta total')
    % colorbar
    % caxis([0 max(power_theta_ncr./power_theta_tot)])
    % 
    % subplot(2,4,3)
    % topoplot(power_alpha_ncr./power_alpha_tot,'19can.locs');
    % title('Potencia alpha total')
    % colorbar
    % caxis([0 max(power_alpha_ncr./power_alpha_tot)])
    % 
    % subplot(2,4,4)
    % topoplot(power_beta_ncr./power_beta_tot,'19can.locs');
    % title('Potencia beta total')
    % colorbar
    % caxis([0 max(power_beta_ncr./power_beta_tot)])
    % 
    % subplot(2,4,5)
    % topoplot(power_gammabaja_ncr./power_gammabaja_tot,'19can.locs');
    % title('Potencia Gamma Baja')
    % colorbar
    % caxis([0 max(power_gammabaja_ncr./power_gammabaja_tot)])
    % 
    % subplot(2,4,6)
    % topoplot(power_gammaalta_ncr./power_gammaalta_tot,'19can.locs');
    % title('Potencia Gamma Alta')
    % colorbar
    % caxis([0 max(power_gammaalta_ncr./power_gammaalta_tot)])
    % 
    % subplot(2,4,7)
    % topoplot(power_ripple1_ncr./power_ripple1_tot,'19can.locs');
    % title('Potencia Ripple1')
    % colorbar
    % caxis([0 max(power_ripple1_ncr./power_ripple1_tot)])
    % 
    % subplot(2,4,8)
    % topoplot(power_ripple2_ncr./power_ripple2_tot,'19can.locs');
    % title('Potencia Ripple 2')
    % colorbar
    % caxis([0 max(power_ripple2_ncr./power_ripple2_tot)])
    % 
    % saveas(fig4,[results_dir filesep 'fig_nocrisis2.png'])
end
