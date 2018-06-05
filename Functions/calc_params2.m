function params = calc_params2(Pxx,frecs,Fs,frec_band)



p_band=Pxx(frecs>frec_band(1) & frecs<frec_band(2));
        
params(1)=medfreq(Pxx,Fs,[frec_band(1) frec_band(2)]); %F_median
params(2)=rms(p_band); % RMS
params(3)=SampEn(2 ,0.2*std(p_band),p_band); % Sample Entropy


