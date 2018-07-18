% Create vector:

subject='PAT_7-EEG_1050';

results_path='/Users/carolinamigliorelli/Documents/MATLAB/TFG-SEIZURE/results';

feature_names_f={'Power','delta','theta','alpha','beta','gammalow','gammahigh','ripple1','ripple2','entropy','rms','mfrec'};


load([results_path subject filesep 'features.mat'])

features_promedio=permute(mean(features,2),[3 1 2]);


features_all=[features_promedio; es_crisis_win];feature_names_f={'Power','delta','theta','alpha','beta','gammalow','gammahigh','ripple1','ripple2','entropy','rms','mfrec'};

features_all=features_all';