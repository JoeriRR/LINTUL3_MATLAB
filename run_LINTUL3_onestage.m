%% add subfolder and clear variables
addpath('Init_scripts')
clear variables; close all; clc;
%% global parameters
SetGlobal;
%% get weather and season data
year_of_data = 1987;
station = 6;
yod = num2str(year_of_data);
weather_file =['ICWEATHR\NLD', num2str(station),'.', yod(2:4)];
%% season params
global DOYEM
% the day of the year on which crop emerges
DOYEM = 90;
% start time simulation
STTIME = 1;
% max season length
season_max_length = 280;
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
%% initial conditions
[SLA,FLVT,NBALAN,TRANRF,TAGBM,RNSOIL,NTAC,NTAG,PEVAP,WATBAL,RWRT,RWST,LUECAL,RWSO,NDEMTO,RWLVG, FRTWET,FSTT,FSOT,FRACT,TSUM,NRF,TEVAP,TTRAN,TRAIN,TRUNOF,TIRRIG,TDRAIN,TEXPLO,CUMPAR,GTSUM,WDRT,CBALAN,NLOSSR,NLOSSL,FERTNS,FERTN,TNSOIL,NUPTT,LAI,WLV,WLVG,WLVD,WSO,WST,WRT,ROOTD,WA,WC,WCCR,DVS,NNI,ANLV,ANST,ANRT,ANSO,NUPTR]  =  deal(zeros(241,1));
SetInit;
%% Simulate
k = 1;
TTSUM = TSUMAN+TSUMMT;
while k<season_max_length && TSUM(k) < TTSUM && DVS(k) < 2.01
    [TSUM(k+1),ROOTD(k+1),WA(k+1),WC(k+1),WCCR(k+1),TEXPLO(k+1),TEVAP(k+1),TTRAN(k+1),TRUNOF(k+1),TIRRIG(k+1),TDRAIN(k+1),TRAIN(k+1),DVS(k+1),NNI(k+1),SLA(k+1),LAI(k+1),NDEMTO(k+1),TNSOIL(k+1),NUPTT(k+1),ANLV(k+1),ANST(k+1),ANRT(k+1),ANSO(k+1),NLOSSL(k+1),NLOSSR(k+1),WLVG(k+1),WLVD(k+1),WST(k+1),WSO(k+1),WRT(k+1),TAGBM(k+1),NTAC(k+1),LUECAL(k+1),CUMPAR(k+1), GTSUM(k+1),WDRT(k+1),DSLR,DVS1,DVS2]...
        = LINTUL3_onestage(DOY(k),LAI(k),TSUM(k),DTEFF(k),WC(k),ROOTD(k),RAIN(k),DAVTMP(k),VP(k),DTR(k),WN(k),TEXPLO(k),WA(k),TEVAP(k),TTRAN(k),TRUNOF(k),TIRRIG(k),TDRAIN(k),TRAIN(k),WLVG(k),WST(k),ANLV(k),ANST(k),ANRT(k),ANSO(k),WRT(k),WSO(k),TNSOIL(k),NUPTT(k),NLOSSL(k),NLOSSR(k),WLVD(k),WDRT(k),GTSUM(k),CUMPAR(k),DSLR,DVS1,DVS2,0,1,0,1);
    k = k+1;
end
%% compare
[time,~,DVS_fst,TSUM_fst,TAGBM_fst,WST_fst,WLVG_fst,WLVD_fst,WSO_fst,LAI_fst,NTAC_fst,WRT_fst,GTSUM_fst, CBALAN_fst, TRANRF_fst, NNI_fst, SLA_fst, FRACT_fst,FRTWET_fst,FLVT_fst,FSTT_fst,FSOT_fst,RWLVG_fst,RWST_fst,RWRT_fst,RWSO_fst,CUMPAR_fst,LUECAL_fst,NUPTT_fst,TTRAN_fst,TEVAP_fst,PEVAP_fst,NBALAN_fst,WATBAL_fst,NUPTR_fst,TNSOIL_fst,NDEMTO_fst,RNSOIL_fst,FERTN_fst,FERTNS_fst,WA_fst,TIRRIG_fst,TRAIN_fst,TEXPLO_fst,TRUNOF_fst,TDRAIN_fst] = importfile('res.dat'); 
figure(1)
plot(1:k-1,LAI(1:k-1))
hold on
plot(1:k-1,LAI_fst(1:k-1))
xlabel('day k')
legend('LAI one stage', 'LAI fst')

MSE_LAI = norm(LAI(1:240)-LAI_fst,2)/length(1:240);
MSE_TEVAP = norm(TEVAP(1:240)-TEVAP_fst,2)/length(1:240);










