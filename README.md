# Muscle-Fatigue-Experiment
Experiment uses EMG signals to detect muscle fatigue.
Above is the MATLAB file used to analyze EMG signals of the biceps to figure out when the biceps fatigue occurs.  
To figur out if muscle fatigue occured, the EMG signal will show that there is a decrease in frequency and increase
in amplitude. 
The signal was cleaned by applying a bandpass filter with a cutoff frequency of 5 Hz a high pass filter and 350 Hz with 
a low pass filter using the butterworth's filter in MATLAB. The sampling frequency of 1000Hz was used to have a
Nyquist rate og 2000Hz. Moreover, the noise from the powerline of 60 Hz was also removed by using a 
band stop filter with a cutoff frequency of 58 to 62 Hz. After applying a fitler, the signal was rectified 
by taking the absolute value and subtracting it from the mean of the filtered 
EMG data. Then, an envolope of the EMG signal was created  using the rms() function. 
After the enveloe, an FFT was taken to determine if the mean and median frequency of the FFT before the 
fatigue and after fatigue shifted. 


