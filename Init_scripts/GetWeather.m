%% season params
global DOYEM
% the day of the year on which crop emerges
DOYEM = 90;
% start time simulation
STTIME = 88;
% max season length
season_max_length = 360;
season_start_date(1:2) = monthday(STTIME); % season start: [M,D], e.g.: [3,4] equals March 4th
DOY = STTIME:season_max_length-1;          % day of year

%% weather data
[DTR,RAIN,TN,TX,~,VP,WN] = Read_Weatherfile(weather_file);
days = DOY(STTIME):DOY(end);
DTR = DTR(days); RAIN = RAIN(days); TN = TN(days); TX = TX(days); VP = VP(days); WN = WN(days);
% DTR is daily total irradiation in MJ/m^2
% RAIN=max(0,RAIN+dRAIN(1:242)');
DAVTMP  = (TN+TX)/2;           % daily average temperature
DTEFF   = max(0,DAVTMP-TBASE); % daily effective temperature