% Create vector:

patient={'PAT_7','PAT_7','PAT_9','PAT_9','PAT_12','PAT_12','PAT_13','PAT_13'};
filenames={'EEG_1059','EEG_1050','EEG_1275','EEG_1220','EEG_1692','EEG_1681','EEG_1933','EEG_1923'};


results_path=[pwd filesep 'results'];
feature_names_f={'Power','delta','theta','alpha','beta','gammalow','gammahigh','ripple1','ripple2','entropy','rms','mfrec'};

features_4clases=[];
features_2clases=[];
features_30sec=[];
for i = 1:length(patient) % de esta manera se calculan caracteristicas para TODOS los pacientes
    % para hacer un detector personalizado para cada paciente:
    % for i = 1:2 --> PAT_7
    % for i = 3:4 --> PAT_9
    % for i = 5:6 --> PAT_12
    % for i = 7:8 --> PAT_13
        
       load([results_path patient{i} '-' filenames{i} filesep 'features2.mat'])
       features_promedio=permute(mean(features,2),[3 1 2]);
       features_std=permute(std(permute(features,[2 1 3])),[3 2 1]);
       aux=[features_promedio;features_std; es_crisis_win];
        
       es_crisis_win(es_crisis_win~=1)=0;
        
       aux2=[features_promedio; features_std; es_crisis_win];
       features_4clases=[features_4clases aux];
       features_2clases=[features_2clases aux2];
        
        
       load([results_path patient{i} '-' filenames{i} filesep 'features_30s.mat'])
      
       features_promedio=permute(mean(Matrix_mean_segmento,2),[3 1 2]);
       features_std=permute(mean(Matrix_std_segmento,2),[3 1 2]);
       aux=[features_promedio;features_std; label_segmento'];
        
       features_30sec=[features_30sec aux];
   
end

features_4clases=features_4clases';
features_2clases=features_2clases';
features_30sec=features_30sec';

