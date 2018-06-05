function params = calc_params(Pxx,frecs,Fs,frec_band)



p_band=Pxx(frecs>frec_band(1) & frecs<frec_band(2));
        
params(1)=frecs(Pxx==max(p_band)); % F_peak
params(2)=medfreq(Pxx,Fs,[frec_band(1) frec_band(2)]); %F_median
params(3)=var(p_band); % Varianza
params(4)=rms(p_band); % RMS
params(5)=SampEn(2 ,0.2*std(p_band),p_band); % Sample Entropy
params(6)=kurtosis(p_band); % Kurtosis
params(7)=skewness(p_band); % Skewness


