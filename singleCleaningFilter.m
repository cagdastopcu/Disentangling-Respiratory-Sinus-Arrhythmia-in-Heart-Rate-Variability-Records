clear all;
close all;
clc,

%%
%%Author: Cagdas Topcu
%% 2018 Potsdam and Weiz
%% Reference
% Topçu, Ç., Frühwirth, M., Moser, M., Rosenblum, M., & Pikovsky, A. (2018). Disentangling Respiratory Sinus Arrhythmia in Heart %%%%Rate Variability Records. arXiv preprint arXiv:1802.00683.
format long g

pathPhaseDataset = 'D:\Projekte\2419_MUG COSMOS\Daten\Phase_data\';

tic
subjecNumber = 25;%%7 added again

varDiff = zeros(1,subjecNumber);
var_ksi_dot = zeros(1,subjecNumber);
var_phi_e_dot = zeros(1,subjecNumber);
var_Qcalc = zeros(1,subjecNumber);
SDNN = zeros(subjecNumber,3);
RMSSD = zeros(subjecNumber,3);
meanRR = zeros(subjecNumber,3);

NN50 = zeros(subjecNumber,3);
pNN50 = zeros(subjecNumber,3);
vagalTone = zeros(subjecNumber,3);
varianceIntervals = zeros(subjecNumber,3);
intervalsAll = {};

for subjectNo=1:subjecNumber;
%     iStr = mat2str(i);%10-46 arasý
    iStr = {'01' '02' '03' '04' '05' '06' '08' '09' '10' ...
        '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' ...
        '21' '22' '23' '24' '25' '26'};
%     iStr = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' ...
%         '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' ...
%         '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' ...
%         '31' '32' '33' '34' '35' '36' '37' '38' '39' '40' ...
%         '41' '42' '43' '44' '45' '46'};    


name_e = ['phi_e_' iStr{subjectNo} '.mat'];
name_r = ['phi_r_' iStr{subjectNo} '.mat'];
name_Q = ['Q_' iStr{subjectNo} '.mat'];
load(fullfile(pathPhaseDataset, name_e))
load(fullfile(pathPhaseDataset, name_r))
load(fullfile(pathPhaseDataset, name_Q))

Q_PhaseECG = Q(:,1);
Q_PhaseResp = Q(:,2);
Q_submeanQ = Q(:,4);
Q = Q(:,3);
phi_e = wrapTo2Pi(phi_e);
phi_r = wrapTo2Pi(phi_r);

phi_e_old = phi_e;
phi_r_old = phi_r;

%% first three derivatives of phi_e
dt = 0.001;
n_phi_e = length(phi_e);
t = (dt.*(0:n_phi_e))';
fsample = 1000;
norder=4;   % order of the fitting polynomial
sl=4 ;      % window semi-length
wl=2*sl+1;  % window length
[b,g] = sgolay(norder,wl);   % Calculate S-G coefficients
phi_e=unwrap(phi_e);
phi_e=phi_e(:); phi_e_dot=phi_e; 
for n=sl+1:length(phi_e)-sl
    phi_e_dot(n)=g(:,2)'*phi_e(n-sl:n+sl);
end
phi_e_dot=phi_e_dot(sl+1:end-sl)*fsample;
phi_e=phi_e(sl+1:end-sl);   % Truncating the phase in order to
phi_e_old=phi_e_old(sl+1:end-sl);
% synchronize it with the derivative               
%% first three derivatives of phi_r
dt = 0.001;
n_phi_r = length(phi_r);
t = (dt.*(0:n_phi_r))';
fsample = 1000;
norder=4;   % order of the fitting polynomial
sl=4 ;      % window semi-length
wl=2*sl+1;  % window length
[b,g] = sgolay(norder,wl);   % Calculate S-G coefficients
phi_r=unwrap(phi_r);
phi_r=phi_r(:); phi_r_dot=phi_r; 
for n=sl+1:length(phi_r)-sl
    phi_r_dot(n)=g(:,2)'*phi_r(n-sl:n+sl);
end
phi_r_dot=phi_r_dot(sl+1:end-sl)*fsample;
phi_r=phi_r(sl+1:end-sl);   % Truncating the phase in order to
phi_r_old=phi_r_old(sl+1:end-sl);

Qcalc = zeros(1,length(phi_e_old));
F = scatteredInterpolant(Q_PhaseECG,Q_PhaseResp,Q_submeanQ);
for i = 1:length(phi_e_old)
% Qcalc(i) = griddata(QPhaseECG,QPhaseResp,Q,phi_e_old(i),phi_r_old(i));
Qcalc(i) = F(phi_e_old(i),phi_r_old(i));
end
% Qcalc(isnan(Qcalc))=0;
% interp2(q1,phi_e(10),phi_r(10))
ksi_dot = phi_e_dot-(Qcalc');

% % figure,
% % 
% % plot(phi_e_dot)
% % hold on
% % plot(ksi_dot,'r')
% % plot(Qcalc,'k')
% % I = legend('$\dot{\varphi}$','$\dot{\xi}$','Q','fontsize',18)
% % set(I,'interpreter','latex');
% % set(I,'FontSize',18);
% % xlabel('Sample (N)','fontsize',18)
% % ylabel('rad/s','fontsize',18)
%% -0.0107587184048687
% varienceDiff = var(ksi_dot-mean(ksi_dot))-(var(phi_e_dot-mean(phi_e_dot))-var(Qcalc-mean(Qcalc)))
var_ksi_dot(subjectNo) = var(ksi_dot);
var_phi_e_dot(subjectNo) = var(phi_e_dot);
var_Qcalc(subjectNo) = var(Qcalc);
varDiff(subjectNo) = var(ksi_dot)-(var(phi_e_dot)-var(Qcalc));

% dim = [.2 .5 .3 .3];
% varienceDiffStr = num2str(varienceDiff);
% str = ['Comparison of ' '$\dot{\varphi}$, ' '$\dot{\xi}$, and ' 'Q'];
% % I1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
% % set(I,'interpreter','latex');
% tSI = title(str,'fontsize',20)
% set(tSI,'interpreter','latex');
currentFolder = pwd;
mkdir(currentFolder,'results\firsttask')
fileNameFirstTaskvarDiff = [currentFolder '\' 'results' '\' 'firsttask' '\'...
    'varDiff.mat'];
fileNameFirstTaskvar_phi_e_dot = [currentFolder '\' 'results' '\' 'firsttask' '\'...
    'var_phi_e_dot.mat'];
fileNameFirstTaskvar_ksi_dot = [currentFolder '\' 'results' '\' 'firsttask' '\'...
    'var_ksi_dot.mat'];
fileNameFirstTaskvar_Qcalc = [currentFolder '\' 'results' '\' 'firsttask' '\'...
    'var_Qcalc.mat'];
save(fileNameFirstTaskvarDiff,'varDiff')
save(fileNameFirstTaskvar_phi_e_dot,'var_phi_e_dot')
save(fileNameFirstTaskvar_ksi_dot,'var_ksi_dot')
save(fileNameFirstTaskvar_Qcalc,'var_Qcalc')

ksi = cumsum(ksi_dot)/1000;
% ksi = cumtrapz(ksi_dot,[0.001:0.001:(0.001*length(ksi_dot))]);

phi_e = unwrap(phi_e_old);
phi_r = unwrap(phi_r_old);
% plot(phi_e)
% hold on
% plot(phi_r,'r')
%% phi_e has 306 selected cycle floor(max(phi_e)/(2*pi))
v_phi_e = 1:length(phi_e);
selected_cycle_number = floor(max(phi_e)/(2*pi)); %306
piVectorphi_e = 2:2:selected_cycle_number*2;
piVectorphi_e = pi.*piVectorphi_e;
vq_phi_e = interp1(phi_e,v_phi_e,piVectorphi_e,'spline');

%% ksi probably similar

v_phi_e = 1:length(phi_e);

piVectorphi_e = 2:2:selected_cycle_number*2;
piVectorphi_e = pi.*piVectorphi_e;
vq_ksi = interp1(ksi,v_phi_e,piVectorphi_e,'spline');

phi_eVariability = diff(vq_phi_e);
%last RR interval was removed
phi_eVariability = phi_eVariability(1:end-1);
ksiVariability = diff(vq_ksi);
ksiVariability = ksiVariability(1:end-1);
%% Q probably similar
 
QcalcNorm = Qcalc + (mean(phi_e_dot)-mean(Qcalc));
% QcalcNorm = Qcalc;

QcalcNorm = cumsum(QcalcNorm)/1000;
v_Qcalc = 1:length(QcalcNorm);

piVectorQcalc = 2:2:selected_cycle_number*2;
piVectorQcalc = pi.*piVectorQcalc;
% QcalcNorm = Qcalc + (mean(phi_e)-mean(Qcalc));
vq_Qcalc = interp1(QcalcNorm,v_Qcalc,piVectorQcalc,'spline');

QcalcVariability = diff(vq_Qcalc);
QcalcVariability = QcalcVariability(1:end-1);
% QcalcVariability = QcalcVariability + (mean(phi_eVariability)-mean(QcalcVariability));
%% Respiration intervals

%% Respiration 129 cycle

% v_phi_r = 1:length(phi_r);
% 
% piVectorphi_r = 2:2:129*2;
% piVectorphi_r = pi.*piVectorphi_r;
% vq_phi_r = interp1(phi_r,v_phi_r,piVectorphi_r,'spline');
% 
% phi_rVariability = diff(vq_phi_r);

%% plotting
% figure,
% plot(phi_eVariability,'r.')
% 
% hold on
% 
% plot(ksiVariability,'b*')
% % plot(phi_rVariability,'ko')
% plot(QcalcVariability,'ko')
% I = legend('Intervals of $\varphi_e$','Intervals of $\xi$','Intervals of Q');
% set(I,'interpreter','latex');
% set(I,'FontSize',18);
% 
% xlabel('Peaks','fontsize',18)
% ylabel('msec','fontsize',18)

phi_eIntervals = phi_eVariability;
ksiIntervals = ksiVariability;
QIntervals = QcalcVariability;
%% classical HRV
%% Variance
variance_phi_e = var(phi_eIntervals); %msec
variance_ksi= var(ksiIntervals); %msec
variance_Q = var(QIntervals); %msec

varianceIntervals(subjectNo,:) = [variance_phi_e variance_ksi variance_Q];
% varianceIntervals.aa = [variance_phi_e variance_ksi variance_Q];

%% Standard deviation of RR

SDNN_phi_e = std(phi_eIntervals); %msec

SDNN_ksi = std(ksiIntervals); %msec

SDNN_Q = std(QIntervals); %msec

SDNN(subjectNo,:) = [SDNN_phi_e SDNN_ksi SDNN_Q];
%% Vagal Trend

vagalTrend_phi_e = diff(phi_eIntervals); %msec.

vagalTrend_ksi = diff(ksiIntervals); %msec.

vagalTrend_Q = diff(QIntervals); %msec.
%% Vagal Tone

vagalTone_phi_e = log10(abs(median(vagalTrend_phi_e)));

vagalTone_ksi = log10(abs(median(vagalTrend_ksi)));

vagalTone_Q = log10(abs(median(vagalTrend_Q)));

vagalTone(subjectNo,:) = [vagalTone_phi_e vagalTone_ksi vagalTone_Q];

%% mean RR and mean HR

meanRR_phi_e = mean(phi_eIntervals); %msec.
meanHR_phi_e = mean(60000./phi_eIntervals); %b/min

meanRR_ksi = mean(ksiIntervals); %msec.
meanHR_ksi = mean(60000./ksiIntervals); %b/min

meanRR_Q = mean(QIntervals); %msec.
meanHR_Q = mean(60000./QIntervals); %b/min

meanRR(subjectNo,:) = [meanRR_phi_e meanRR_ksi meanRR_Q];

%% RMSSD

RMSSD_phi_e = sqrt(mean(diff(phi_eIntervals).*diff(phi_eIntervals)));

RMSSD_ksi = sqrt(mean(diff(ksiIntervals).*diff(ksiIntervals)));

RMSSD_Q = sqrt(mean(diff(QIntervals).*diff(QIntervals)));

RMSSD(subjectNo,:) = [RMSSD_phi_e RMSSD_ksi RMSSD_Q];

%% NN50
alpha = 50; %ms

NN50_phi_e = sum( abs(diff(phi_eIntervals)) >= alpha);

NN50_ksi = sum( abs(diff(ksiIntervals)) >= alpha);

NN50_Q = sum( abs(diff(QIntervals)) >= alpha);


NN50(subjectNo,:) = [NN50_phi_e NN50_ksi NN50_Q];

%% PPNAlpha PPN50
alpha = 50; %ms

pNN50_phi_e = sum( abs(diff(phi_eIntervals)) >= alpha)/length(diff(phi_eIntervals));

pNN50_ksi = sum( abs(diff(ksiIntervals)) >= alpha )/length(diff(ksiIntervals));

pNN50_Q = sum( abs(diff(QIntervals)) >= alpha )/length(diff(QIntervals));

pNN50(subjectNo,:) = [pNN50_phi_e pNN50_ksi pNN50_Q];

%% logRSA

logRSA_phi_e = log(median(abs(diff(phi_eIntervals))));

logRSA_ksi = log(median(abs(diff(ksiIntervals))));

logRSA_Q = log(median(abs(diff(QIntervals))));

logRSA(subjectNo,:) = [logRSA_phi_e logRSA_ksi logRSA_Q];

%% files
mkdir(currentFolder,'results\classicalHRV')
fileNameClassicalHRVSDNN = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'SDNN.mat'];
fileNameClassicalHRVvagalTone = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'vagalTone.mat'];
fileNameClassicalHRVmeanRR = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'meanRR.mat'];
fileNameClassicalHRVRMSSD = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'RMSSD.mat'];
fileNameClassicalHRVNN50 = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'NN50.mat'];
fileNameClassicalHRVlogRSA = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'logRSA.mat'];
fileNameClassicalHRVpNN50 = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'pNN50.mat'];
fileNameClassicalHRVvarianceIntervals = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'varianceIntervals.mat'];

save(fileNameClassicalHRVSDNN,'SDNN')
save(fileNameClassicalHRVvagalTone,'vagalTone')
save(fileNameClassicalHRVmeanRR,'meanRR')
save(fileNameClassicalHRVRMSSD,'RMSSD')
save(fileNameClassicalHRVNN50,'NN50')
save(fileNameClassicalHRVpNN50,'pNN50')
save(fileNameClassicalHRVlogRSA,'logRSA')
save(fileNameClassicalHRVvarianceIntervals,'varianceIntervals')
% % % % % % % % % end
% % % % % % % % % toc
%% nonlinear analysis - entropy measures
rate_phi_e = 60000./phi_eIntervals;
rate_ksi = 60000./ksiIntervals;
rate_Q = 60000./QIntervals;

winsizeRate = length(rate_phi_e);
winincRate = winsizeRate;
datawinRate = hanning(winsizeRate);
dispstatus = 0;

ApEn_rate_phi_e(subjectNo) = getapenfeat(rate_phi_e',winsizeRate,winincRate,datawinRate,dispstatus);
ApEn_rate_ksi(subjectNo) = getapenfeat(rate_ksi',winsizeRate,winincRate,datawinRate,dispstatus);
ApEn_rate_Q(subjectNo) = getapenfeat(rate_Q',winsizeRate,winincRate,datawinRate,dispstatus);


SampEn_rate_phi_e(subjectNo) = getsampenfeat(rate_phi_e',winsizeRate,winincRate,datawinRate,dispstatus);
SampEn_rate_ksi(subjectNo) = getsampenfeat(rate_ksi',winsizeRate,winincRate,datawinRate,dispstatus);
SampEn_rate_Q(subjectNo) = getsampenfeat(rate_Q',winsizeRate,winincRate,datawinRate,dispstatus);


%% frequency analysis
fs = 7;
Ts = 1000./fs;

rate_phi_e = 60000./phi_eIntervals;
% rate_phi_r = 60000./phi_rIntervals;
rate_ksi = 60000./ksiIntervals;
rate_Q = 60000./QIntervals;
%% X:\COSMOS\Cagdas\Matlab-tbx\hrv
% now 7 Hz
phi_e_rate_4Hz = regsampling(rate_phi_e,phi_eIntervals, Ts, 'pchip');
% phi_r_rate_4Hz = regsampling(rate_phi_r,phi_rIntervals, 250, 'pchip');
ksi_rate_4Hz = regsampling(rate_ksi,ksiIntervals, Ts, 'pchip');
Q_rate_4Hz = regsampling(rate_Q,ksiIntervals, Ts, 'pchip');

%Settings: Window length 57.14 s
% window = 400;
% noverlap = window/4;
% f = [0:0.01:0.5];
%Settings: Window length 300 s
expected_length = 410*fs;
window = 300*fs;
noverlap = 2*window - expected_length;
if length(phi_e_rate_4Hz) < expected_length
	if length(phi_e_rate_4Hz) < window
		error('%s: segment length too short: %d samples', iStr{subjectNo}, length(phi_e_rate_4Hz));
	end	
	noverlap = 2*window - (length(phi_e_rate_4Hz) - 14);	% set to actual length minus 2 seconds
	warning('%s: noverlap reduced from %g to %g', iStr(subjectNo), 2*window - expected_length, noverlap);
end
f = window;

% [ps_phi_e_rate_4Hz, freq_phi_e_rate_4Hz] = getFFT4Hz(phi_e_rate_4Hz',winsize,wininc,datawin,dispstatus);
[ps_phi_e_rate_4Hz, freq_phi_e_rate_4Hz] = pwelch(phi_e_rate_4Hz-mean(phi_e_rate_4Hz),window,noverlap,f,fs);
% [ps_phi_r_rate_4Hz, freq_phi_r_rate_4Hz] = pwelch(phi_r_rate_4Hz,window,noverlap,f,fs);
% [ps_ksi_rate_4Hz, freq_ksi_rate_4Hz] = getFFT4Hz(ksi_rate_4Hz',winsize,wininc,datawin,dispstatus);
[ps_ksi_rate_4Hz, freq_ksi_rate_4Hz] = pwelch(ksi_rate_4Hz-mean(ksi_rate_4Hz),window,noverlap,f,fs);
[ps_Q_rate_4Hz, freq_Q_rate_4Hz] = pwelch(Q_rate_4Hz-mean(Q_rate_4Hz),window,noverlap,f,fs);

% [ps_Q_rate_4Hz, freq_Q_rate_4Hz] = getFFT4Hz(Q_rate_4Hz',winsize,wininc,datawin,dispstatus);
% [ps_Q_rate_4Hz, freq_Q_rate_4Hz] = pwelch(Q_rate_4Hz,window,noverlap,f,fs);
[VLF_phi_e(subjectNo), LF_phi_e(subjectNo), HF_phi_e(subjectNo), LFHFratio_phi_e(subjectNo),nLF_phi_e(subjectNo), nHF_phi_e(subjectNo)] = frequencyFeatures(freq_phi_e_rate_4Hz,ps_phi_e_rate_4Hz); 
[VLF_ksi(subjectNo), LF_ksi(subjectNo), HF_ksi(subjectNo), LFHFratio_ksi(subjectNo), nLF_ksi(subjectNo), nHF_ksi(subjectNo)] = frequencyFeatures(freq_ksi_rate_4Hz,ps_ksi_rate_4Hz); 
[VLF_Q(subjectNo), LF_Q(subjectNo), HF_Q(subjectNo), LFHFratio_Q(subjectNo), nLF_Q(subjectNo), nHF_Q(subjectNo)] = frequencyFeatures(freq_Q_rate_4Hz,ps_Q_rate_4Hz); 


ps_phi_e_rate_4Hz_All(subjectNo,:) = ps_phi_e_rate_4Hz;
ps_ksi_rate_4Hz_All(subjectNo,:) = ps_ksi_rate_4Hz;
ps_Q_rate_4Hz_All(subjectNo,:) = ps_Q_rate_4Hz;
ps_phi_ksi_Q_All = [ps_phi_e_rate_4Hz_All ps_ksi_rate_4Hz_All ps_Q_rate_4Hz_All];

% phi_eIntervals = phi_eVariability;
% ksiIntervals = ksiVariability;
% QIntervals = QcalcVariability;
intervalsAll{subjectNo} = [phi_eIntervals; ksiIntervals; QIntervals];

% rates4HzAll(subjectNo,:) = [phi_e_rate_4Hz ksi_rate_4Hz Q_rate_4Hz];

clear data
clear phi_r
clear phi_e
clear Q
end
%% nonl
ApEn = [ApEn_rate_phi_e' ApEn_rate_ksi' ApEn_rate_Q'];
SampEn = [SampEn_rate_phi_e' SampEn_rate_ksi' SampEn_rate_Q'];
VLF = [VLF_phi_e' VLF_ksi' VLF_Q'];
LF = [LF_phi_e' LF_ksi' LF_Q'];
HF = [HF_phi_e' HF_ksi' HF_Q'];
LFHFratio = [LFHFratio_phi_e' LFHFratio_ksi' LFHFratio_Q'];
nLF = [nLF_phi_e' nLF_ksi' nLF_Q'];
nHF = [nHF_phi_e' nHF_ksi' nHF_Q'];

fileNameClassicalHRVVLF = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'VLF.mat'];
fileNameClassicalHRVLF = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'LF.mat'];
fileNameClassicalHRVHF = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'HF.mat'];
fileNameClassicalHRVLFHFratio = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'LFHFratio.mat'];
fileNameClassicalHRVnLF = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'nLF.mat'];
fileNameClassicalHRVnHF = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'nHF.mat'];
fileNameClassicalHRVApEn = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'ApEn.mat'];
fileNameClassicalHRVSampEn = [currentFolder '\' 'results' '\' 'classicalHRV' '\'...
    'SampEn.mat'];

save(fileNameClassicalHRVVLF,'VLF')
save(fileNameClassicalHRVLF,'LF')
save(fileNameClassicalHRVHF,'HF')
save(fileNameClassicalHRVLFHFratio,'LFHFratio')
save(fileNameClassicalHRVnLF,'nLF')
save(fileNameClassicalHRVnHF,'nHF')

save(fileNameClassicalHRVApEn,'ApEn')
save(fileNameClassicalHRVSampEn,'SampEn')

mkdir(currentFolder,'results\PowerSpectrum')
fileNameps = [currentFolder '\' 'results' '\' 'PowerSpectrum' '\'...
    'ps_phi_ksi_Q_All.mat'];
save(fileNameps,'ps_phi_ksi_Q_All')
% fileNamefreq = mkdir(currentFolder,'results\PowerSpectrum')
fileNamefreq = [currentFolder '\' 'results' '\' 'PowerSpectrum' '\'...
    'freq.mat'];
freq = freq_phi_e_rate_4Hz;
save(fileNamefreq,'freq')
mkdir(currentFolder,'results\Intervals')
fileNameintervalsAll = [currentFolder '\' 'results' '\' 'Intervals' '\'...
    'intervalsAll.mat'];
save(fileNameintervalsAll,'intervalsAll')
% fileNamerates4HzAll = [currentFolder '\' 'results' '\' 'Intervals' '\'...
%     'rates4HzAll.mat'];
% save(fileNamerates4HzAll,'rates4HzAll')

% plot(freq_phi_r_rate_4Hz,10*log10(ps_phi_r_rate_4Hz),'m','linewidth',2);

toc
