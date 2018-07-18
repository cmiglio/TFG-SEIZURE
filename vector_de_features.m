% Create vector:

subject='PAT_7-EEG_1050';

results_path='C:\Users\Caro\Documents\MATLAB\TFGs\TFG-SEIZURE\results';

feature_names_f={'Power','delta','theta','alpha','beta','gammalow','gammahigh','ripple1','ripple2','entropy','rms','mfrec'};


load([results_path filesep subject])

features_promedio=permute(mean(features,2),[3 1 2]);


features_all=[features_promedio; es_crisis_win];feature_names_f={'Power','delta','theta','alpha','beta','gammalow','gammahigh','ripple1','ripple2','entropy','rms','mfrec'};

