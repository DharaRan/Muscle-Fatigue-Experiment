%Developer: Dhara Rana
%Contact: djr32@njit.edu
%Developed on:Oct 18 2016
%Version 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% STEP 1 Load DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
load('Trial2_Exp.mat')
Trial= 'Trial 2';
emg=data-0.5;% EMG data
fs=1000; % sampling frequency
t=1/fs:1/fs:length(data)/fs;
figure(1)
plot(t,emg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% %%%%%%% STEP 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Filter EMG Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% pass a band pass filter by%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Low pass filter: fc of 500 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%high pass filter: fc of 5hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%High Pass filter
fcH= 5;
WnH= (fcH)/(fs/2);

[bH,aH]=butter(4,WnH,'high');
emgHigh=filtfilt(bH,aH,emg);

%Filter out 60 Hz noise from power
fc60_1=58;
fc60_2=62;
Wn1= (fc60_1)/(fs/2);
Wn2=(fc60_2)/(fs/2);
[b60,a60]=butter(4,[Wn1 Wn2],'stop');
emgband=filtfilt(b60,a60,emgHigh);
%Filter out 120 Hz noise from power
fc60_1=118;
fc60_2=122;
Wn1= (fc60_1)/(fs/2);
Wn2=(fc60_2)/(fs/2);
[b60,a60]=butter(4,[Wn1 Wn2],'stop');
emgband=filtfilt(b60,a60,emgband);
%Filter out 180 Hz noise from power
fc60_1=178;
fc60_2=182;
Wn1= (fc60_1)/(fs/2);
Wn2=(fc60_2)/(fs/2);
[b60,a60]=butter(4,[Wn1 Wn2],'stop');
emgband=filtfilt(b60,a60,emgband);
%Filter out 240 Hz noise from power
fc60_1=238;
fc60_2=242;
Wn1= (fc60_1)/(fs/2);
Wn2=(fc60_2)/(fs/2);
[b60,a60]=butter(4,[Wn1 Wn2],'stop');
emgband=filtfilt(b60,a60,emgband);
%Filter out 300 Hz noise from power
fc60_1=298;
fc60_2=302;
Wn1= (fc60_1)/(fs/2);
Wn2=(fc60_2)/(fs/2);
[b60,a60]=butter(4,[Wn1 Wn2],'stop');
emgband=filtfilt(b60,a60,emgband);
%Filter out 360 Hz noise from power
fc60_1=358;
fc60_2=362;
Wn1= (fc60_1)/(fs/2);
Wn2=(fc60_2)/(fs/2);
[b60,a60]=butter(4,[Wn1 Wn2],'stop');
emgband=filtfilt(b60,a60,emgband);

%Filter out 420 Hz noise from power
fc60_1=418;
fc60_2=422;
Wn1= (fc60_1)/(fs/2);
Wn2=(fc60_2)/(fs/2);
[b60,a60]=butter(4,[Wn1 Wn2],'stop');
emgband=filtfilt(b60,a60,emgband);

%Low pass filter
fcL= 450;
Wn= (fcL)/(fs/2);
[b,a]=butter(4,Wn,'low');
filtemg=filtfilt(b,a,emgband);

hold on
plot(t,filtemg,'g')

xlabel('Time (sec)')
ylabel('Potential (mV)')
title(strcat(Trial,':Filtered EMG Signal'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% STEP 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Split up EMG signal into parts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EMG before the onset of muscle contraction (0-5 sec)%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deg0 = 7.09;%time in seconds recorded when no muscle contracts
endofRest=deg0/(1/fs);
trial1_0 = filtemg(1:endofRest);%it starts at to end 7.08 9 sec
trial1_length0 = length(trial1_0);%length when it reaches that trial
trial1_time_0 = (1:1:endofRest)./fs;

%EMG start of muscle contraction to time of 30 degrees(5sec to t of 30 degree sec)
deg30 = 95;%time in seconds recorded when muscle starts to fatigue
Atdeg0=endofRest+1;
beginFatgue=deg30/(1/fs);
trial1_30 = filtemg(Atdeg0:beginFatgue);%it starts at 70
trial1_length30 = length(trial1_30);%length when it reaches that trial
trial1_time_30 = (Atdeg0:1:beginFatgue)./fs;

%EMG start of muscle fatigue(40 degrees) to time of 70 degrees(t of 30 degree sec to 70 degree)
deg70 = 158;%time recorded when muscle starts to fatigue
Atfatigue30=beginFatgue+1;
enddeg70=deg70/(1/fs);
trial1_70 = filtemg(Atfatigue30:enddeg70);%it starts at 95
trial1_length70 = length(trial1_70);%length when it reaches that trial
trial1_time_70 = (Atfatigue30:1:enddeg70)./fs;

%Time of 70 degrees to 5 seconds after(t of 70 degree to last 5 seconds)
backNormalPos = deg70+deg0; %time recorded when muscle fatigue ends
At70=enddeg70+1;
endnormalpos=backNormalPos/(1/fs);
trial1_Norm = filtemg(At70:endnormalpos);%it starts at 95
trial1_lengthNorm = length(trial1_Norm);%length when it reaches that trial
trial1_time_Norm = (At70:1:endnormalpos)./fs;

figure(123)
p1=plot(trial1_time_0,trial1_0,'Color', [0.5 0 0.5]);
hold on
p2=plot(trial1_time_30,trial1_30,'Color',[0 0.5 0]);
p3=plot(trial1_time_70,trial1_70,'Color',[0.5 0 0]);
p4=plot(trial1_time_Norm,trial1_Norm,'b');

xlabel('Time (sec)')
ylabel('Potential (mV)')
title(strcat(Trial,':Filtered EMG Signal'))
legend([p1 p2 p3 p4],{'Before Muscle Contraction','During Muscle Contraction','During Muscle Fatigue','After Muscle Fatigue'},'Location','southeast')
ylim([-0.1 0.15]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% STEP 4: Rectify and RMS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rectify the data and then Create linear envelope of signal
%Linear envelope using a simple Root Mean Square Moving Average Filter - see Matlabs exmaple of moving average filters in the help file
%Try different window sizes to get a smooth envelope that clearly gives the
%shape of the data
%retifying the raw emg data for each section
rEMG=abs(filtemg);
EMG_line1 = abs(trial1_0); 
EMG_line2 = abs(trial1_30);
EMG_line3 = abs(trial1_70);
EMG_line4 = abs(trial1_Norm);
figure(224)
p11= plot(trial1_time_0,EMG_line1,'Color', [0.5 0 0.5]);
hold on
p12= plot(trial1_time_30,EMG_line2,'Color',[0 0.5 0]);
p13= plot(trial1_time_70,EMG_line3,'Color',[0.5 0 0]);
p14= plot(trial1_time_Norm,EMG_line4,'b');

rmsenvelop1=rms384(rEMG(1:endnormalpos),3000,2999,1);

meanAMP_BC =mean(rmsenvelop1(1:endofRest));
meanAMP_OC_BF =mean(rmsenvelop1(Atdeg0:beginFatgue));
meanAMP_DF =mean(rmsenvelop1(Atfatigue30:enddeg70));
meanAMP_AF =mean(rmsenvelop1(At70:endnormalpos));

meanAMP =[meanAMP_BC meanAMP_OC_BF meanAMP_DF meanAMP_AF]

pRMS1=plot(t(1:endnormalpos),rmsenvelop1,'k');

xlabel('Time (sec)')
ylabel('Potential (mV)')
title(strcat(Trial,':Filtered Rectified EMG Signal'))
legend([p11 p12 p13 p14 pRMS1],{'Before Muscle Contraction','Muscle Contraction','Muscle Fatigue','Rest','Linear Envelope'},'Location','northeast')
ylim([0 0.15]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% STEP 5 FFT and MFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%FFT 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%FFT before Contraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(200)
BCN = length(filtemg(1:endofRest));
BCindex = 0:1:BCN-1;%Range of frequency (starting at DC freq)
BCfres = fs/BCN; 
BCfreq = BCindex.*BCfres; %Determine frequency resolution
BCfiltemg_FFT= fft(filtemg(1:endofRest));%fft of force
BCfiltemg_FFT = abs(BCfiltemg_FFT);% absolute value of fft
BCfiltemg_FFT = abs(BCfiltemg_FFT/(length(BCfiltemg_FFT)/2));% absolute value of fft
BCfiltemg_FFT = BCfiltemg_FFT(1:length(BCfiltemg_FFT)/2);
BCfres = fs/BCN; 
BCindex = 0:1:(BCN-1)/2;
BCfreq = BCindex.*BCfres;

p200=plot(BCfreq,BCfiltemg_FFT,'Color', [0.5 0 0.5]);
xlabel('Frequency (Hz)')
ylabel('Potential (mV)')
title(strcat(Trial,':Frequency Spectrum Before Muscle Contraction'))
ylim([0 0.006])

% Power Spectrum
[pxx, f] = periodogram(filtemg(1:endofRest), [], [], fs);
BCMPF=medfreq(pxx,f)
meanBCMPF=meanfreq(pxx,f)
hold on
figure(202)
plot(f, 10 * log10(pxx),'Color', [0.5 0 0.5])
l201=line([BCMPF BCMPF],[-160 20],'LineWidth',2.5,'Color',[1 0.5 0]);
L202=line([meanBCMPF meanBCMPF],[-160 20],'LineWidth',2.5,'Color','k');
legend([l201,L202],{strcat('Median Power Frequency (', num2str(BCMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meanBCMPF),' Hz)')})
xlabel('Frequency (Hz)')
ylabel('Power(mV^2/Hz)')
title(strcat(Trial,':Power Spectrum Before Muscle Contraction'))

figure(200)
hold on
l201=line([BCMPF BCMPF],[0 0.006],'LineWidth',2.5,'Color',[1 0.5 0]);
L202=line([meanBCMPF meanBCMPF],[0 0.006],'LineWidth',2.5,'Color','k');
legend([p200, l201,L202],{'FFT of Before Contraction',strcat('Median Power Frequency (', num2str(BCMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meanBCMPF),' Hz)')})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%FFT Onset contraction (OC) and before Fatigue(BF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(300)
N = length(filtemg(Atdeg0:beginFatgue));
index = 0:1:N-1;%Range of frequency (starting at DC freq)
fres = fs/N; 
freq = index.*fres; %Determine frequency resolution
filtemg_FFT= fft(filtemg(Atdeg0:beginFatgue));%fft of force
OCBFfiltemg_FFT = abs(filtemg_FFT/(length(filtemg_FFT)/2));% absolute value of fft
OCBFfiltemg_FFT = OCBFfiltemg_FFT(1:length(OCBFfiltemg_FFT)/2);
fres = fs/N; 
index = 0:1:(N-1)/2;
freq = index.*fres;
p300=plot(freq,OCBFfiltemg_FFT,'Color',[0 0.5 0]);
xlabel('Frequency (Hz)')
ylabel('Potential (mV)')
title(strcat(Trial,':Frequency Spectrum During Muscle Contraction'))
ylim([0 0.0015])

% Power Spectrum
[pxx, f] = periodogram(filtemg(Atdeg0:beginFatgue), [], [], fs);
figure(302)
plot(f, 10 * log10(pxx),'Color',[0 0.5 0])
xlabel('Frequency (Hz)')
ylabel('Power(mV^2/Hz)')
title(strcat(Trial,':Power Spectrum During Muscle Contraction'))
OCBFMPF=medfreq(pxx,f)
meanOCBFMPF=meanfreq(pxx,f)
hold on
l302=line([OCBFMPF OCBFMPF],[-160 20],'LineWidth',1.5,'Color',[1 .5 0]);
L303=line([meanOCBFMPF meanOCBFMPF],[-160 20],'LineWidth',1.5,'Color','k');
legend([l302,L303],{strcat('Median Power Frequency (', num2str(OCBFMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meanOCBFMPF),' Hz)')})

figure(300)
hold on
l302=line([OCBFMPF OCBFMPF],[0 0.0015],'LineWidth',1.5,'Color',[1 .5 0]);
L303=line([meanOCBFMPF meanOCBFMPF],[0 0.0015],'LineWidth',1.5,'Color','k');
legend([p300, l302,L303],{'FFT of During Muscle Contraction',strcat('Median Power Frequency (', num2str(OCBFMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meanOCBFMPF),' Hz)')})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %%FFT During Fatigue(BF)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(400)
DFN = length(filtemg(Atfatigue30:enddeg70));
DFindex = 0:1:DFN-1;%Range of frequency (starting at DC freq)
dffres = fs/DFN; 
dffreq = DFindex.*dffres; %Determine frequency resolution
DFfiltemg_FFT= fft(filtemg(Atfatigue30:enddeg70));%fft of force
DFfiltemg_FFT = abs(DFfiltemg_FFT);% absolute value of fft
DFfiltemg_FFT = abs(DFfiltemg_FFT/(length(DFfiltemg_FFT)/2));% absolute value of fft
DFfiltemg_FFT = DFfiltemg_FFT(1:length(DFfiltemg_FFT)/2);
dffres = fs/DFN; 
DFindex = 0:1:(DFN-1)/2;
dffreq = DFindex.*dffres;

p400=plot(dffreq,DFfiltemg_FFT,'Color',[0.5 0 0]);
xlabel('Frequency (Hz)')
ylabel('Potential (mV)')
title(strcat(Trial,':Frequency Spectrum During Muscle Fatigue'))
ylim([0 0.002])


% Power Spectrum
[pxx, f] = periodogram(filtemg(Atfatigue30:enddeg70), [], [], fs);
figure(402)
plot(f, 10 * log10(pxx),'Color',[0.5 0 0])
xlabel('Frequency (Hz)')
ylabel('Power(mV^2/Hz)')
title(strcat(Trial,':Power Spectrum During Muscle Fatigue'))
dfMPF=medfreq(pxx,f)
meandfMPF=meanfreq(pxx,f)
hold on
L402=line([dfMPF dfMPF],[-160 20],'LineWidth',1.5,'Color',[1 .5 0]);
L403=line([meandfMPF meandfMPF],[-160 20],'LineWidth',1.5,'Color','k');
legend([L402,L403],{strcat('Median Power Frequency (', num2str(dfMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meandfMPF),' Hz)')})

figure(400)
hold on
L402=line([dfMPF dfMPF],[0 0.002],'LineWidth',1.5,'Color',[1 .5 0]);
L403=line([meandfMPF meandfMPF],[0 0.002],'LineWidth',1.5,'Color','k');
legend([p400, L402,L403],{'FFT of During Muscle Fatigue',strcat('Median Power Frequency (', num2str(dfMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meandfMPF),' Hz)')})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %%FFT After Muscle Fatigue(BF)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(500)
AFN = length(filtemg(At70:endnormalpos));
AFindex = 0:1:AFN-1;%Range of frequency (starting at DC freq)
AFfres = fs/AFN; 
AFfreq = AFindex.*AFfres; %Determine frequency resolution
AFfiltemg_FFT= fft(filtemg(At70:endnormalpos));%fft of force
AFfiltemg_FFT = abs(AFfiltemg_FFT);% absolute value of fft
AFfiltemg_FFT = abs(AFfiltemg_FFT/(length(AFfiltemg_FFT)/2));% absolute value of fft
AFfiltemg_FFT = AFfiltemg_FFT(1:length(AFfiltemg_FFT)/2);
AFfres = fs/AFN; 
AFindex = 0:1:(AFN-1)/2;
AFfreq = AFindex.*AFfres;

p500=plot(AFfreq,AFfiltemg_FFT,'b');
xlabel('Frequency (Hz)')
ylabel('Potential (mV)')
title(strcat(Trial,':Frequency Spectrum: After Muscle Fatigue'))
ylim([0 0.0015])

% Power Spectrum
[pxx, f] = periodogram(filtemg(At70:endnormalpos), [], [], fs);
figure(502)
plot(f, 10 * log10(pxx),'b')
xlabel('Frequency (Hz)')
ylabel('Power(mV^2/Hz)')
title(strcat(Trial,':Power Spectrum: After Muscle Fatigue'))
AFMPF=medfreq(pxx,f)
meanAFMPF=meanfreq(pxx,f)
hold on
L502=line([AFMPF AFMPF],[-140 20],'LineWidth',1.5,'Color',[1 .5 0]);
L503=line([meanAFMPF meanAFMPF],[-140 20],'LineWidth',1.5,'Color','k');
legend([L502,L503],{strcat('Median Power Frequency (', num2str(AFMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meanAFMPF),' Hz)')})

figure(500)
hold on
L502=line([AFMPF AFMPF],[0 0.0015],'LineWidth',1.5,'Color',[1 .5 0]);
L503=line([meanAFMPF meanAFMPF],[0 0.0015],'LineWidth',1.5,'Color','k');
legend([p500, L502,L503],{'FFT of After Muscle Fatigue',strcat('Median Power Frequency (', num2str(AFMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meanAFMPF),' Hz)')})

MPF=[BCMPF OCBFMPF dfMPF AFMPF]
MeanPF=[meanBCMPF meanOCBFMPF meandfMPF meanAFMPF]


%Ploting FFT in subplots
figure(90)

subplot(4,1,1)

    p200=plot(BCfreq,BCfiltemg_FFT,'Color', [0.5 0 0.5]);
    hold on
    l201=line([BCMPF BCMPF],[0 0.006],'LineWidth',2.5,'Color',[1 .5 0]);
    L202=line([meanBCMPF meanBCMPF],[0 0.006],'LineWidth',2.5,'Color','k');
    legend([p200, l201,L202],{'FFT of Before Contraction',strcat('Median Power Frequency (', num2str(BCMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meanBCMPF),' Hz)')})
    xlabel('Frequency (Hz)')
    ylabel('Potential (mV)')
    title(strcat(Trial,':Frequency Spectrum Before Muscle Contraction'))
    ylim([0 0.006])

subplot(4,1,2)
    p300=plot(freq,OCBFfiltemg_FFT,'Color',[0 0.5 0]);
    hold on
    l302=line([OCBFMPF OCBFMPF],[0 0.0015],'LineWidth',1.5,'Color',[1 .5 0]);
    L303=line([meanOCBFMPF meanOCBFMPF],[0 0.0015],'LineWidth',1.5,'Color','k');
    legend([p300, l302,L303],{'FFT of During Muscle Contraction',strcat('Median Power Frequency (', num2str(OCBFMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meanOCBFMPF),' Hz)')})
    xlabel('Frequency (Hz)')
    ylabel('Potential (mV)')
    title(strcat(Trial,':Frequency Spectrum During Muscle Contraction'))
    ylim([0 0.0015])

subplot(4,1,3)
    p400=plot(dffreq,DFfiltemg_FFT,'r');
    hold on
    L402=line([dfMPF dfMPF],[0 0.003],'LineWidth',1.5,'Color',[1 .5 0]);
    L403=line([meandfMPF meandfMPF],[0 0.003],'LineWidth',1.5,'Color','k');
    legend([p400, L402,L403],{'FFT of During Muscle Fatigue',strcat('Median Power Frequency (', num2str(dfMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meandfMPF),' Hz)')})
    xlabel('Frequency (Hz)')
    ylabel('Potential (mV)')
    title(strcat(Trial,':Frequency Spectrum During Muscle Fatigue'))
    ylim([0 0.003])
    
subplot(4,1,4)
    p500=plot(AFfreq,AFfiltemg_FFT,'b');
    hold on
    L502=line([AFMPF AFMPF],[0 0.0015],'LineWidth',1.5,'Color','[1 .5 0]');
    L503=line([meanAFMPF meanAFMPF],[0 0.0015],'LineWidth',1.5,'Color','k');
    legend([p500, L502,L503],{'FFT of After Muscle Fatigue',strcat('Median Power Frequency (', num2str(AFMPF),' Hz)'),strcat('Mean Power Frequency (', num2str(meanAFMPF),' Hz)')})
    xlabel('Frequency (Hz)')
    ylabel('Potential (mV)')
    title(strcat(Trial,':Frequency Spectrum After Muscle Fatigue'))
    ylim([0 0.0015])

