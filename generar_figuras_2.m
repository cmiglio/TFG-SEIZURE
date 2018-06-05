clear all
close all
addpath('functions');
patient={'PAT_7','PAT_7','PAT_9','PAT_9','PAT_12','PAT_12','PAT_13','PAT_13'};
filenames={'EEG_1059','EEG_1050','EEG_1275','EEG_1220','EEG_1692','EEG_1681','EEG_1933','EEG_1923'};
feature_names={'Power','\delta','\theta','\alpha','\beta','{\gamma}_{low}','{\gamma}_{high}','ripple_{1}','ripple_{2}','entropy','rms','mfrec'};
feature_names_f={'Power','delta','theta','alpha','beta','gammalow','gammahigh','ripple1','ripple2','entropy','rms','mfrec'};


cm=colormap(hot);
cm(1:10,:)=[];
cm=cm(size(cm,1):-1:1,:);
cm2=colormap(redblue);
mkdir('results_figures')


for i=1:length(patient)
    load(['results' patient{i} '-' filenames{i} filesep 'features.mat'])
    filepos='19can.locs';    
    if strcmp(patient{i},'PAT_9')
        features(:,[5 13],:)=[]; % eliminar canales artefactuados
        filepos='17can.locs';
    elseif strcmp(patient{i},'PAT_7')
        features(:,14,:)=[]; % eliminar canales artefactuados
        filepos='18can.locs';
    end
    
    n_crisis=sum(es_crisis_win==1);
    inter_index=find(es_crisis_win==0);
    i_idx=inter_index(randi(length(inter_index),n_crisis,1));
    
    for j=1:size(features,3)
        % Crisis
        feat_crisis(:,j)=mean(features(es_crisis_win==1,:,j));
        % Baseline
        feat_bas(:,j)=mean(features(:,:,j));
        % Interictal
        feat_inter(:,j)=mean(features(i_idx,:,j));
        % Pre ictal
        feat_pre(:,j)=mean(features(es_crisis_win==2,:,j));   
        % Post ictal
        feat_pos(:,j)=mean(features(es_crisis_win==3,:,j));   
    end
    
    % Topografias power, entropy rms mfreq
    for j=1:size(features,3)
        
        cmax=max([feat_inter(:,j); feat_pre(:,j); feat_crisis(:,j); feat_pos(:,j)]);
        cmin=min([feat_inter(:,j); feat_pre(:,j); feat_crisis(:,j); feat_pos(:,j)]);
        
        if j==1
            feat_inter(:,j)=10*log10(feat_inter(:,j));
        end
        
        figure(j)
        subplot(2,4,1)
        topoplot(feat_inter(:,j),filepos);
        title(['Interictal: ' feature_names{j} ' absolute'])
        colorbar
        caxis([cmin cmax])
        colormap(gca,cm)
    
        subplot(2,4,2)
        topoplot(feat_pre(:,j),filepos);
        title(['Preictal: ' feature_names{j} ' absolute'])
        colorbar
        caxis([cmin cmax])
        colormap(gca,cm)
        
        subplot(2,4,3)
        topoplot(feat_crisis(:,j),filepos);
        title(['Ictal: ' feature_names{j} ' absolute'])
        colorbar
        caxis([cmin cmax])
        colormap(gca,cm)
        
        subplot(2,4,4)
        topoplot(feat_pos(:,j),filepos);
        title(['Postictal: ' feature_names{j} ' absolute'])
        colorbar
        caxis([cmin cmax])
        colormap(gca,cm)

        cmax=max([feat_inter(:,j)./feat_bas(:,j); feat_pre(:,j)./feat_bas(:,j); feat_crisis(:,j)./feat_bas(:,j); feat_pos(:,j)./feat_bas(:,j)]);
        cmin=min([feat_inter(:,j)./feat_bas(:,j); feat_pre(:,j)./feat_bas(:,j); feat_crisis(:,j)./feat_bas(:,j); feat_pos(:,j)./feat_bas(:,j)]);
        
        cmax=max(abs(cmax),abs(cmin));
        cmin=-cmax+2;
        
        subplot(2,4,5)
        topoplot(feat_inter(:,j)./feat_bas(:,j),filepos);
        title(['Interictal: ' feature_names{j} ' relative'])
        colorbar
        caxis([cmin cmax])
        colormap(gca,cm2)
    
        subplot(2,4,6)
        topoplot(feat_pre(:,j)./feat_bas(:,j),filepos);
        title(['Preictal: ' feature_names{j} ' relative'])
        colorbar
        caxis([cmin cmax])
        colormap(gca,cm2)
        
        subplot(2,4,7)
        topoplot(feat_crisis(:,j)./feat_bas(:,j),filepos);
        title(['Ictal: ' feature_names{j} ' relative'])
        colorbar
        caxis([cmin cmax])
        colormap(gca,cm2)
        
        subplot(2,4,8)
        topoplot(feat_pos(:,j)./feat_bas(:,j),filepos);
        title(['Postictal: ' feature_names{j} ' relative'])
        colorbar
        caxis([cmin cmax])
        colormap(gca,cm2)
    end
    
    
    for j=1:size(features,3)
        fig=figure(j);
        fig.Color=[1 1 1];
        fig.Position=[335 426 1504 564];
        
        saveas(fig,['results_figures' filesep  patient{i} '-' filenames{i} '-' feature_names_f{j} '.fig'])
        saveas(fig,['results_figures' filesep  patient{i} '-' filenames{i} '-' feature_names_f{j} '.png'])
    end  
    
    close all
    clear feat_crisis feat_pos feat_pre feat_inter feat_bas
end