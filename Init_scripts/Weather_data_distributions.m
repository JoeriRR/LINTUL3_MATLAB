%close all
clc
clear variables
%% Verified version
%% Season parameters
global DOYEM
DOYEM = 90;
STTIME = 1;
season_start_date(1:2) = monthday(STTIME); % season start: [M,D], e.g.: [3,4] equals March 4th
season_max_length   = 300;
DOY = STTIME:season_max_length; % day of year
savedata = false;
%% Crop parameters
TBASE   = 0.0;   % [degC] base temperature
%% Weather data
year_of_data = [1966:1998]; %available range of weather files
station = 6;
% init = 46 + days(datetime(year_of_data,season_start_date(1),season_start_date(2))-datetime('1-1-1951'));
% [YYYYMMDD,FG,TG,TN,TX,Q0,UG,UX,UN,EV24] = ...
%     ParseWeatherFile('etmgeg_370.txt',init,init+season_max_length-1);
for i = 1:length(year_of_data)
    yod = num2str(year_of_data(i));    
    weather_file =['C:\Users\s168210\Desktop\Stage\Stage_Joeri\LINTUL\ICWEATHR\NLD', num2str(station),'.', yod(2:4)];
    if station ==1
        [DTR,RAIN,TN,TX,~,VP,WN] = LINTUL_weatherfile_all(weather_file,DOY(1),DOY(end));
    elseif station == 6
        [DTR,RAIN,TN,TX,~,VP,WN] = Read_Weatherfile(weather_file);
        DTRm(:,i) = DTR(DOY(1):DOY(end));
        RAINm(:,i) = RAIN(DOY(1):DOY(end));
        TNm(:,i) = TN(DOY(1):DOY(end));
        TXm(:,i) = TX(DOY(1):DOY(end));
        VPm(:,i) = VP(DOY(1):DOY(end));
        WNm(:,i) = WN(DOY(1):DOY(end));
    end
    % DTR is daily total irradiation in MJ/m^2
    
    % RAIN=max(0,RAIN+dRAIN(1:242)');
    TRAINm(:,i)   = cumsum(RAINm(:,i));        % cumulative rain
    DAVTMPm(:,i)  = (TNm(:,i)+TXm(:,i))/2;     % daily average temperature
    DTEFFm(:,i)   = max(0,DAVTMPm(:,i)-TBASE); % daily effective temperature
end

%% calculate mean and std
mean_DTR = mean(DTRm,2);
mean_RAIN = mean(RAINm,2);
mean_TN = mean(TNm,2);
mean_TX = mean(TXm,2);
mean_VP = mean(VPm,2);
mean_WN = mean(WNm,2);
mean_TRAIN = mean(TRAINm,2);
mean_DAVTMP = mean(DAVTMPm,2);
mean_DTEFF = mean(DTEFFm,2);

std_DTR = std(DTRm,0,2);
std_RAIN = std(RAINm,0,2);
std_TN = std(TNm,0,2);
std_TX = std(TXm,0,2);
std_VP = std(VPm,0,2);
std_WN = std(WNm,0,2);
std_TRAIN = std(TRAINm,0,2);
std_DAVTMP = std(DAVTMPm,0,2);
std_DTEFF = std(DTEFFm,0,2);
%% implement forecast errors
if(savedata == true)
    save('Expected_weather.mat','std_DTR','std_RAIN','std_DAVTMP','std_DTEFF','std_VP','std_WN', 'mean_DTR','mean_RAIN','mean_DAVTMP','mean_DTEFF','mean_VP','mean_WN');
end
%% calculate correlation.
CorrMat = [mean_DTR mean_RAIN mean_DAVTMP mean_VP mean_WN];
Correlation = zeros(size(CorrMat,2),size(CorrMat,2));
for i = 1:size(CorrMat,2)
   for j = 1:size(CorrMat,2)
       R = corrcoef(CorrMat(:,i), CorrMat(:,j));
       Correlation(i,j) = R(1,2);       
   end
end
%% figures
figure(1)
subplot(511)
plot(DOY,std_DTR)
hold on
plot(DOY,mean_DTR)
xlim([DOY(1) DOY(end)])
legend('STD','MEAN')
ylabel('DTR')
title(['Standard deviation STD and MEAN of weatherdata years: ' num2str(year_of_data(1)) '- ' num2str(year_of_data(end))])

subplot(512)
plot(DOY,std_RAIN)
hold on
plot(DOY,mean_RAIN)
xlim([DOY(1) DOY(end)])
ylabel('RAIN')

subplot(513)
plot(DOY,std_VP)
hold on
plot(DOY,mean_VP)
xlim([DOY(1) DOY(end)])
ylabel('VP')

subplot(514)
plot(DOY,std_DAVTMP)
hold on
plot(DOY,mean_DAVTMP)
xlim([DOY(1) DOY(end)])
ylabel('DAVTMP')

subplot(515)
plot(DOY,std_WN)
hold on
plot(DOY,mean_WN)
xlim([DOY(1) DOY(end)])
ylabel('WN')
xlabel('Day of year')

%% example weather data usage during prediction
figure(2)
n_w = 7;
t = 20;
rain = zeros(season_max_length,1);
rain(t:t+n_w) = DTRm(t:t+n_w,4);
rain(t+n_w+1:end) = mean_DTR(t+n_w+1:end);
plot(t:season_max_length,rain(t:season_max_length),t:season_max_length,DTRm(t:season_max_length,4),'Linewidth',2);
xlim([20 100])
xlabel('t [days]')
ylabel('DTR')
title(['weather data iteration step t = ' num2str(t)]);
line([t+n_w t+n_w],[0 25],'color',[0 0 0],'Linestyle','- -');
ylim([0 25])
legend('weatherdata with average', 'real weather data year 1970','forecast horizon');





