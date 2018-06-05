

patient={'PAT_7','PAT_7','PAT_9','PAT_9','PAT_12','PAT_12','PAT_13','PAT_13'};
filenames={'EEG_1059','EEG_1050','EEG_1275','EEG_1220','EEG_1692','EEG_1681','EEG_1933','EEG_1923'};
feature_names={'Power','\delta','\theta','\alpha','\beta','{\gamma}_{low}','{\gamma}_{high}','ripple_{1}','ripple_{2}','entropy','rms','mfrec'};

fig_names={'seizure_abs_feat','seizure_rel_feat','interictal_abs_feat','interictal_rel_feat','pre_abs_feat','pre_rel_feat','pos_abs_feat','pos_rel_feat',...
    'seizure_abs_bands','seizure_rel_bands','interictal_abs_bands','interictal_rel_bands','pre_abs_bands','pre_rel_bands','pos_abs_bands','pos_rel_bands'};

for i=1:length(patient)
    load(['results' patient{i} '-' filenames{i} filesep 'features.mat'])
    
% Valores medios durante la crisis

    for j=1:size(features,3)
        % Crisis
        feat_crisis(:,j)=mean(features(es_crisis_win==1,:,j));
        % Baseline
        feat_bas(:,j)=mean(features(:,:,j));
        % Interictal
        feat_inter(:,j)=mean(features(es_crisis_win==0,:,j));
        % Pre ictal
        feat_pre(:,j)=mean(features(es_crisis_win==2,:,j));   
        % Post ictal
        feat_pos(:,j)=mean(features(es_crisis_win==3,:,j));   
    end

    % Topografias power, entropy rms mfreq
    
    i_feat=[1 10 11 12]; % indices
    for j=1:length(i_feat)
        figure(1)
        subplot(2,2,j)
        topoplot(feat_crisis(:,i_feat(j)),'19can.locs');
        title(['Seizure: ' feature_names{i_feat(j)} ' absolute'])
        colorbar
        caxis([0 max(feat_crisis(:,i_feat(j)))])
    
        figure(2)
        subplot(2,2,j)
        topoplot(feat_crisis(:,i_feat(j))./feat_bas(:,i_feat(j)),'19can.locs');
        title(['Seizure: ' feature_names{i_feat(j)} ' relative'])
        colorbar
        caxis([0 max(feat_crisis(:,i_feat(j))./feat_bas(:,i_feat(j)))])
                
        figure(3)
        subplot(2,2,j)
        topoplot(feat_inter(:,i_feat(j)),'19can.locs');
        title(['Interictal: ' feature_names{i_feat(j)} ' absolute'])
        colorbar
        caxis([0 max(feat_inter(:,i_feat(j)))])
    
        figure(4)
        subplot(2,2,j)
        topoplot(feat_inter(:,i_feat(j))./feat_bas(:,i_feat(j)),'19can.locs');
        title(['Interictal: ' feature_names{i_feat(j)} ' relative'])
        colorbar
        caxis([0 max(feat_inter(:,i_feat(j))./feat_bas(:,i_feat(j)))]) 
    
        figure(5)
        subplot(2,2,j)
        topoplot(feat_pre(:,i_feat(j)),'19can.locs');
        title(['Preictal: ' feature_names{i_feat(j)} ' absolute'])
        colorbar
        caxis([0 max(feat_pre(:,i_feat(j)))])
    
        figure(6)
        subplot(2,2,j)
        topoplot(feat_pre(:,i_feat(j))./feat_bas(:,i_feat(j)),'19can.locs');
        title(['Preictal: ' feature_names{i_feat(j)} ' relative'])
        colorbar
        caxis([0 max(feat_pre(:,i_feat(j))./feat_bas(:,i_feat(j)))])
        
        figure(7)
        subplot(2,2,j)
        topoplot(feat_pos(:,i_feat(j)),'19can.locs');
        title(['Postictal: ' feature_names{i_feat(j)} ' absolute'])
        colorbar
        caxis([0 max(feat_pos(:,i_feat(j)))])
    
        figure(8)
        subplot(2,2,j)
        topoplot(feat_pos(:,i_feat(j))./feat_bas(:,i_feat(j)),'19can.locs');
        title(['Postictal: ' feature_names{i_feat(j)} ' relative'])
        colorbar
        caxis([0 max(feat_pos(:,i_feat(j))./feat_bas(:,i_feat(j)))])
    end
    
    i_feat=2:9;
    for j=1:length(i_feat)
        
        figure(9)
        subplot(2,4,j)
        topoplot(feat_crisis(:,i_feat(j)),'19can.locs');
        title(['Seizure: ' feature_names{i_feat(j)} ' absolute'])
        colorbar
        caxis([0 max(feat_crisis(:,i_feat(j)))])
    
        figure(10)
        subplot(2,4,j)
        topoplot(feat_crisis(:,i_feat(j))./feat_bas(:,i_feat(j)),'19can.locs');
        title(['Seizure: ' feature_names{i_feat(j)} ' relative'])
        colorbar
        caxis([0 max(feat_crisis(:,i_feat(j))./feat_bas(:,i_feat(j)))])
                
        figure(11)
        subplot(2,4,j)
        topoplot(feat_inter(:,i_feat(j)),'19can.locs');
        title(['Interictal: ' feature_names{i_feat(j)} ' absolute'])
        colorbar
        caxis([0 max(feat_inter(:,i_feat(j)))])
    
        figure(12)
        subplot(2,4,j)
        topoplot(feat_inter(:,i_feat(j))./feat_bas(:,i_feat(j)),'19can.locs');
        title(['Interictal: ' feature_names{i_feat(j)} ' relative'])
        colorbar
        caxis([0 max(feat_inter(:,i_feat(j))./feat_bas(:,i_feat(j)))]) 
    
        figure(13)
        subplot(2,4,j)
        topoplot(feat_pre(:,i_feat(j)),'19can.locs');
        title(['Preictal: ' feature_names{i_feat(j)} ' absolute'])
        colorbar
        caxis([0 max(feat_pre(:,i_feat(j)))])
    
        figure(14)
        subplot(2,4,j)
        topoplot(feat_pre(:,i_feat(j))./feat_bas(:,i_feat(j)),'19can.locs');
        title(['Preictal: ' feature_names{i_feat(j)} ' relative'])
        colorbar
        caxis([0 max(feat_pre(:,i_feat(j))./feat_bas(:,i_feat(j)))])
        
        figure(15)
        subplot(2,4,j)
        topoplot(feat_pos(:,i_feat(j)),'19can.locs');
        title(['Postictal: ' feature_names{i_feat(j)} ' absolute'])
        colorbar
        caxis([0 max(feat_pos(:,i_feat(j)))])
    
        figure(16)
        subplot(2,4,j)
        topoplot(feat_pos(:,i_feat(j))./feat_bas(:,i_feat(j)),'19can.locs');
        title(['Postictal: ' feature_names{i_feat(j)} ' relative'])
        colorbar
        caxis([0 max(feat_pos(:,i_feat(j))./feat_bas(:,i_feat(j)))])
    end
    
    mkdir(['results' patient{i} '-' filenames{i} filesep 'figures'])   

    for j=1:size(features,3)
        fig=figure(j);
        fig.Color=[1 1 1];
        if j>8
            fig.Position=[178 126 1502 688];
        end
        saveas(fig,['results' patient{i} '-' filenames{i} filesep 'figures' filesep fig_names{j} '.fig'])
        saveas(fig,['results' patient{i} '-' filenames{i} filesep 'figures' filesep fig_names{j} '.png'])
    end  
    
    close all
end