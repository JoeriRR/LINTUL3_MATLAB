% todo :
% - include globals in functions
% - rescale state variables
% - comments
% - one stage simulator

close all; clear all;clc;
%% Season parameters
global DOYEM
DOYEM = 90;
STTIME = 1;
season_max_length = 280;
season_start_date(1:2) = monthday(STTIME); % season start: [M,D], e.g.: [3,4] equals March 4th
DOY = STTIME:season_max_length-1; % day of year

%% global parameters not grouped yet
global SLAC TBASE NMAXSO DELT WMFAC RDRNS DVSNT DVSNLT WCSUBS DVSDR LSNR FRNX LRNR WCI NFRLVI NFRSTI NFRRTI RGRL LAICR TSUMAN TSUMMT TSUMAG RDRSHM LUE K ROOTDM RRDMAX WCAD WCWET WCST WCWP WCFC TRANCO DRATE IRRIGF TCNT RDRRT RNFLV RNFST RNFRT FNTRT NLUE NLAI NSLA NPART
SLAC = 0.022;
TBASE = 0;
RGRL = 0.009;
LAICR = 4;
RDRSHM = 0.03;
LUE = 2.8;
K = 0.6;
ROOTDM = 1.2;
RRDMAX = 0.012;
WCAD = 0.1; WCWP = 0.2; WCFC = 0.4; WCWET = 0.45; WCST = 0.5;
TRANCO = 8;
DRATE = 30;
IRRIGF = 1;
TCNT = 10;
RDRRT = 0.03;
RNFLV = 0.004; RNFST = 0.002; RNFRT = 0.002;
FNTRT = 0.15;
NLUE = 0.2;
NLAI = 1;
NSLA = 1;
NPART = 1;
WCI = 0.4;
TSUMAN = 800; TSUMMT = 1030; TSUMAG = 800;
NFRLVI = 0.06;  NFRSTI = 0.03; NFRRTI = 0.03;
NMAXSO = 0.0165;
LSNR = 0.50; LRNR = 0.50;
FRNX = 0.5;
DVSNLT = 1.0;
DVSDR = 1.0;
WMFAC = 0.0;
WCSUBS= 0.30;
DELT = 1;
RDRNS = 0.03;
DVSNT = 0.8;
%initial conditions
TSUMI = 0;
DVSI = 0;
ROOTDI = 0.1;
WAI = 1000*ROOTDI*WCI;
WLVGI = 2.4; WSTI  = 0.0; WRTLI = 3.6; WSOI  = 0.0;
SLACFI = SLACF(DVSI);
ISLA = SLAC * SLACFI;
LAII = WLVGI * ISLA;
ANLVI = NFRLVI* WLVGI;
ANSTI = NFRSTI* WSTI;
ANRTI = NFRRTI* WRTLI;
ANSOI = 0.;

%% Weather data
year_of_data = 1987;
station = 6;

% init = 46 + days(datetime(year_of_data,season_start_date(1),season_start_date(2))-datetime('1-1-1951'));
% [YYYYMMDD,FG,TG,TN,TX,Q0,UG,UX,UN,EV24] = ...
%     ParseWeatherFile('etmgeg_370.txt',init,init+season_max_length-1);

yod = num2str(year_of_data);
weather_file =['C:\Users\s168210\Documents\MATLAB\stage\GIT\NLD', num2str(station),'.', yod(2:4)];
[DTR,RAIN,TN,TX,~,VP,WN] = Read_Weatherfile(weather_file);
days = DOY(1):DOY(end);
DTR = DTR(days); RAIN = RAIN(days); TN = TN(days); TX = TX(days); VP = VP(days); WN = WN(days);
% DTR is daily total irradiation in MJ/m^2

% RAIN=max(0,RAIN+dRAIN(1:242)');
%TRAIN   = [cumsum(RAIN)];    % cumulative rain
DAVTMP  = (TN+TX)/2;           % daily average temperature
DTEFF   = max(0,DAVTMP-TBASE); % daily effective temperature

%% Initial conditions
% Overview of states
% ============================
% TSUM  = temperature sum [Cd]
% LAI   = leaf area index [-]
% WLV   = biomass of leaves [kg]
% WLVG  = biomass of green leaves [kg]
% WLVD  = biomass of dead leaves [kg]
% WSO   = biomass storage organ [kg]
% WST   = biomass of stem [kg]
% WRT   = biomass roots [kg]
% ROOTD = root depth [m]
% WA    = water in the soil root-zone [mm]
% WC    = water content 

[SLA,FLVT,NBALAN,TRANRF,TAGBM,RNSOIL,NTAC,NTAG,PEVAP,WATBAL,RWRT,RWST,LUECAL,RWSO,NDEMTO,RWLVG, FRTWET,FSTT,FSOT,FRACT,TSUM,NRF,TEVAP,TTRAN,TRAIN,TRUNOF,TIRRIG,TDRAIN,TEXPLO,CUMPAR,GTSUM,WDRT,CBALAN,NLOSSR,NLOSSL,FERTNS,FERTN,TNSOIL,NUPTT,LAI,WLV,WLVG,WLVD,WSO,WST,WRT,ROOTD,WA,WC,WCCR,DVS,NNI,ANLV,ANST,ANRT,ANSO,NUPTR]  =  deal(zeros(241,1));
WLVG(1) = WLVGI;
DVS(1) = DVSI;
ROOTD(1) = ROOTDI;
WA(1) = WAI;
WC(1) = WCI;
WCCR(1) = WCWP+0.01;
LAI(1) = LAII;
WLVG(1) = WLVGI;
WRT(1) = WRTLI;
WST(1) = WSTI;
WSO(1) = WSOI;
ANLV(1) = ANLVI;
ANST(1) = ANSTI;
ANRT(1) = ANRTI;
ANSO(1) = ANSOI;
DSLR = 1;
%% Simulate
k    = 1;
I = zeros(1,season_max_length);
TTSUM = TSUMAN+TSUMMT;
DVS1 = 0; DVS2 = 0;
while k < 241%season_max_length && TSUM(k) < TTSUM && k <= season_max_length
    %% --------- emergence, temperature sum ---------
    EMERG = ((DOY(k)>=(DOYEM)&& WC(k)>WCWP) && -LAI(k)<0 );
    % photosyntatic effect, latitude amersvoort 571.460
    DAYL = ASTRO(DOY(k),571.460);
    if((TSUM(k)-TSUMAN)<0)
        PHOTPF = PHOTTB(DAYL);
    else
        PHOTPF = 1;
    end
    RTSUM = DTEFF(k)*EMERG; 
    RTSUMP = RTSUM * PHOTPF;
    TSUM(k+1) = TSUM(k) + RTSUMP; 
    %% --------- Rooth Depth ---------
    RROOTD = min(RRDMAX*(WC(k) >= WCWP)*EMERG,ROOTDM-ROOTD(k));
    ROOTD(k+1) = ROOTD(k) + RROOTD;
    
    %% --------- Soil moisture balance --------- 
    
    [PEVAP(k),PTRAN] = PENMAN( DAVTMP(k),VP(k),DTR(k),LAI(k),WN(k));    
    [EVAP, TRAN, WCCR(k+1),DSLR] = EVAPTR( PEVAP(k),PTRAN,ROOTD(k),WA(k),WCAD,WCWP,WCFC,WCWET,WCST,TRANCO,DELT,WMFAC,RAIN(k),DSLR);
    
    [DRAIN,RUNOFF,IRRIG] = DRUNIR(RAIN(k),EVAP,TRAN,IRRIGF,DRATE,DELT,WA(k),ROOTD(k),WCFC,WCST,WMFAC);
    % TRANRF verified
    if(PTRAN <= eps)
        TRANRF(k) = TRAN;
    else
        TRANRF(k) = TRAN/PTRAN; 
    end    
    % update water levels
    EXPLOR = 1000*RROOTD*WCSUBS;
    
    RWA = (RAIN(k)+EXPLOR+IRRIG) - (RUNOFF+TRAN+EVAP+DRAIN);  % rainfed
    WA(k+1) = WA(k)+RWA;
    % water content in rootzone
    WC(k+1) = WA(k)*0.001/ROOTD(k);
    % total exploration
    TEXPLO(k+1) = TEXPLO(k) + EXPLOR;
    % crop transipartion and soil evaporation
    TEVAP(k+1) = TEVAP(k) + EVAP;
    TTRAN(k+1) = TTRAN(k) + TRAN;
    % other water balance integrals
    TRUNOF(k+1) = TRUNOF(k) + RUNOFF;
    TIRRIG(k+1) = TIRRIG(k) + IRRIG;
    TDRAIN(k+1) = TDRAIN(k) + DRAIN;
    TRAIN(k+1) = TRAIN(k) + RAIN(k);
    % update water balance
    WATBAL(k) = WA(k)-WAI-TRAIN(k)-TEXPLO(k)-TIRRIG(k)+TRUNOF(k)+TTRAN(k)+TEVAP(k)+TDRAIN(k);
    %% --------- Calculation of Develement stage  DVS and update NNI ---------
    [DVS(k),DVS1,DVS2] = SUBDVS(DOY(k),DOYEM,TSUM(k),TSUMAN,TSUMMT,DVS1,DVS2); % verified

    % maximum N concentration in the leaves, from which the values of the
    % stem and roots are derived, as a function of development stage.
    NMAXLV = NMXLV(DVS(k));
    NMAXST = LSNR * NMAXLV;
    NMAXRT = LRNR * NMAXLV;
    % maximum nitrogen concentration of stem and leaves
    [NOPTLV,NOPTST] = NOPTM(FRNX,NMAXLV,NMAXST);
    % total vegatative biomass
    TBGMR = WLVG(k)+WST(k);
    % maximum N content in the plant
    NOPTS = NOPTST* WST(k);
    NOPTL = NOPTLV* WLVG(k);  
    if(TBGMR == 0)
        NOPTMR = (NOPTL+NOPTS);
    else
        NOPTMR = (NOPTL+NOPTS)/TBGMR;
    end
    NUPGMR = ANLV(k)+ANST(k);
    if(TBGMR == 0)
        NRMR = (WLVG(k)*RNFLV+WST(k)*RNFST);
        NFGMR  = NUPGMR;
    else
        NRMR = (WLVG(k)*RNFLV+WST(k)*RNFST)/(TBGMR);
        NFGMR = NUPGMR/TBGMR;
    end
    NNI(k) = NNINDX(DOY(k),DOYEM,EMERG,NFGMR,NRMR,NOPTMR);
    SLA(k) = SLAC * SLACF(DVS(k))*exp(-NSLA*(1-NNI(k)));
    %% --------- soil nitrogen supply ---------
    RDRTMP = RDRT(DAVTMP(k));
    % growth rate and dry matter production of plant organs
    FRTWET(k) = FRTTB(DVS(k)); % verified
    FLVT(k) = FLVTB(DVS(k)); % verified
    FSTT(k) = FSTTB(DVS(k)); % verified
    FSOT(k) = FSOTB(DVS(k)); % verified
    [FRT,FLV,FST,FSO] = SUBPAR(NPART,TRANRF(k),NNI(k),FRTWET(k),FLVT(k),FSTT(k),FSOT(k));
    FRACT(k) = FLV+FST+FSO+FRT; % verified
    
    %% calculate leaf growth
    [PARINT,GTOTAL] = GROWTH(DOY(k),DOYEM,DTR(k),K,NLUE,LAI(k),LUE,TRANRF(k),NNI(k));
    GLV = FLV * GTOTAL;    
    GLAI = GLA(DOY(k),DOYEM,DTEFF(k),LAII,RGRL,DELT,SLA(k),LAI(k),GLV,NLAI,WC(k),WCWP,DVS(k),TRANRF(k),NNI(k));
    [RDRDV,RDRSH,RDR,DLV,DLVS,DLVNS,DLAIS,DLAINS,DLAI] = DEATHL(DOY(k),DOYEM,TSUM(k),TSUMAG,RDRTMP,RDRSHM,LAI(k),LAICR,WLVG(k),RDRNS,NNI(k),SLA(k));
    RLAI = GLAI - DLAI;        
    LAI(k+1) = LAI(k) + RLAI; % hier gaat t ergens fout
    %% --------- soil nitrogen supply continueed ---------
    
    [DRRT,RNLDLV,RNLDRT] = RNLD(DVS(k),WRT(k),RDRRT,RNFLV,DLV,RNFRT,DVSDR);
    [RWLVG(k),RWRT(k),RWST(k),RWSO(k)] = RELGR(DOY(k),DOYEM,EMERG,WLVGI,WRTLI,WSTI,WSOI,GTOTAL,FLV,FRT,FST,FSO,DLV,DRRT,DELT);
    [NDEML,NDEMS,NDEMR,NDEMSO] = NDEMND(NMAXLV,NMAXST,NMAXRT,NMAXSO,WLVG(k),WST(k),WRT(k),WSO(k),ANLV(k),ANST(k),ANRT(k),ANSO(k),TCNT,DELT);
    NDEMTO(k) = max(0,NDEML+NDEMS+NDEMR); % verified
    
    % Nitrogren limiting factor at low moisture conditions
    if(DVS(k)<DVSNLT && WC(k)>= WCWP)
        NLIMIT = 1;
    else
        NLIMIT = 0;
    end   
    % Soil N supply
    RTMIN  = 0.10 * EMERG * NLIMIT;
    
    % fertilizer application
    NRF(k) = NRFTAB(DOY(k));
    FERTN(k) = FERTAB(DOY(k));
    FERTNS(k) = FERTN(k) * NRF(k); % verified
    % total nutrient uptake
    if(DOY(k)<DOYEM)
        NUPTR(k) = 0;
    else
        NUPTR(k) = max(0., min (NDEMTO(k), TNSOIL(k)))* NLIMIT/ DELT;
    end 
    [RNULV,RNUST,RNURT] = RNUSUB(DOY(k),DOYEM,EMERG,NDEML,NDEMS,NDEMR,NUPTR(k),NDEMTO(k),ANLVI,ANSTI,ANRTI,DELT);
    [ATNLV,ATNST,ATNRT,ATN] = NTRLOC(ANLV(k),ANST(k),ANRT(k),WLVG(k),WST(k),WRT(k),RNFLV,RNFST,RNFRT,FNTRT);
    if(DVS(k)<DVSNT)
        NSUPSO = 0;
    else
        NSUPSO = ATN/TCNT;
    end    
    RNSO = min(NDEMSO,NSUPSO);
    [RNTLV,RNTST,RNTRT] = NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN);
    
    
    %% updating states
    RNLV = RNULV-RNTLV-RNLDLV;
    RNST = RNUST-RNTST; 
    RNRT = RNURT-RNTRT-RNLDRT; 
 
    % Change in organic N in soil
    RNSOIL(k) = FERTNS(k)/DELT -NUPTR(k) + RTMIN;
    % Amount of inorganic N in soil as function of fertizlier
    TNSOIL(k+1) = TNSOIL(k) + RNSOIL(k);
    % Total uptake of N over time (g N m-2)
    NUPTT(k+1) = NUPTT(k) + NUPTR(k);    
    %{
    unused
    % N concentration of the leaves, stem, roots and storage organs
    if(WLVG(k) == 0)
        NFLV = ANLV;
    else
        NFLV = ANLV/WLVG(k);
    end
    if(WST(k) == 0)
        NFST = ANST;
    else
        NFST = ANST/WST(k);
    end
    if(WLVG(k) == 0)
        NFRT = ANRT;
    else
        NFRT = ANRT/WRT(k);
    end
    if(WLVG(k) == 0)
        NFSO = ANSO;
    else
        NFSO = ANSO/WSO(k);
    end    
    %}
     % actual N concent of organs
    ANLV(k+1) = ANLV(k) + RNLV;
    ANST(k+1) = ANST(k) + RNST;
    ANRT(k+1) = ANRT(k) + RNRT;
    ANSO(k+1) = ANSO(k) + RNSO;
    % Nitrogen Balance
    NLOSSL(k+1) = RNLDLV + NLOSSL(k);
    NLOSSR(k+1) = NLOSSR(k) + RNLDRT;
    NBALAN(k) = NUPTT(k)+(ANLVI+ANSTI+ANRTI+ANSOI)-(ANLV(k)+ANST(k)+ANRT(k)+ANSO(k)+NLOSSL(k)+NLOSSR(k)); 
    
    % weight of green leaves, dead leaves, stem, storage orangs and roots
    WLVG(k+1) = WLVG(k)+RWLVG(k);
    WLVD(k+1) = WLVD(k)+DLV;
    WST(k+1) = WST(k)+RWST(k); % verified
    WSO(k+1) = WSO(k)+RWSO(k); 
    WRT(k+1) = WRT(k)+RWRT(k);
    % total leaf weigth
    WLV = WLVG(k)+WLVD(k);
    % total above ground biomass
    TAGBM(k) = WLV+WST(k)+WSO(k); 
    % Water Nitrogren stress factor
    RNW = min(TRANRF(k),NNI(k));
    % Biomass carbon balance
    GTSUM(k+1) = GTSUM(k)+GTOTAL;
    WDRT(k+1) = WDRT(k) + DRRT;
    CBALAN(k) = GTSUM(k)+(WRTLI+WLVGI+WSTI+WSOI)-(WLV+WST(k)+WSO(k)+WRT(k)+WDRT(k));
    
    % For calculation of LUE
    PAR = 0.5*DTR(k);
    LUECAL(k) = GTOTAL/PAR;
    CUMPAR(k+1) = CUMPAR(k) + PAR;
    NTAG(k) = ANLV(k)+ANST(k)+ANSO(k);
    if(NTAG(k) == 0)
        NTAC(k) = NTAG(k);
    else
        NTAC(k) = NTAG(k)/TAGBM(k);
    end
    %% update day
    k = k+1;

end
% resize state variables
%{
DOY     = DOY(1:k);
TSUM    = TSUM(1:k);
LAI     = LAI(1:k);
WLV     = WLV(1:k);
WLVG    = WLVG(1:k);
WLVD    = WLVD(1:k);
WSO     = WSO(1:k);
WST     = WST(1:k);
WRT     = WRT(1:k);
ROOTD   = ROOTD(1:k);
WA      = WA(1:k);
WC      = WC(1:k);
WCCR    = WCCR(1:k);
WSOTHA  = WSO/100;
ANLV    = ANLV(1:k);
%}

save('Z.mat','DVS','TSUM', 'TAGBM','WST', 'WLVG','WLVD','WSO','LAI' ,'NTAC' ,'WRT','GTSUM','CBALAN','TRANRF', 'NNI','SLA', ...
'FRACT','FRTWET','FLVT','FSTT', 'FSOT','RWLVG','RWST','RWRT', 'RWSO','CUMPAR','LUECAL', 'NUPTT','TTRAN','TEVAP','PEVAP','NBALAN', 'WATBAL', ... 
'NUPTR','TNSOIL','NDEMTO','RNSOIL','FERTN','FERTNS','WA','TIRRIG','TRAIN','TEXPLO','TRUNOF','TDRAIN')


%% Comparison to FST
% [TIME_f,LAI_f,WSOTHA_f] =importFSTresult('res.dat');
% hold on
% plot(TIME_f,WSOTHA_f,'--','LineWidth',2);
% plot(TIME_f,LAI_f,'--','LineWidth',2);
% legend('MATLAB: Storage organ biomass [ton/ha]','MATLAB: LAI [-]','FST: Storage organ biomass [ton/ha]','FST: LAI [-]','Location','NW')
%
% WSOTHA_norm = norm(WSOTHA_f'-WSOTHA)
% LAI_norm    = norm(LAI_f'-LAI)
%

%% functions to include photoperiodicity
function [c] = PHOTTB(tau)
PHOTTB_tab = [0,0; 8,1;10,1;12,1;18,1];
c = interp1(PHOTTB_tab(:,1),PHOTTB_tab(:,2),tau);
end
function [c] = RDRT(tau)
RDRT_tab = [-10,0; 10,0.02; 15,0.03; 30,0.05; 50, 0.09];
c = interp1(RDRT_tab(:,1),RDRT_tab(:,2),tau);
end
function [c] = SLACF(tau)
SLACF_tab = [0,1; 2,1; 2.1,1];
c = interp1(SLACF_tab(:,1),SLACF_tab(:,2),tau);
end
function [c] = NMXLV(tau)
NMXLV_tab = [0,0.06;0.4,0.04;0.7,0.03;1,0.02;2,0.014; 2.1,0.014];
c = interp1(NMXLV_tab(:,1),NMXLV_tab(:,2),tau);
end
%% Biomass partitioning
% for leaves (LV), stems (ST), storage organs (SO) and roots (RT) added 2.4
% to make sure no NANS can occur
function [c] = FRTTB(tau)
FRTTB_tab = [0,0.6;0.33,0.58;0.4,0.55;0.8,0.1;1,0;2,0;2.4,0];
c = interp1(FRTTB_tab(:,1),FRTTB_tab(:,2),tau);
end
function [c] = FLVTB(tau)
FLVTB_tab = [ 0, 0.4; 0.33,0.42; 0.4,0.405; 0.8,0.36; 1,0.1; 1.01,0; 2, 0; 2.4,1];
c = interp1(FLVTB_tab(:,1),FLVTB_tab(:,2),tau);
end
function [c] = FSTTB(tau)
FSTTB_tab = [0.,0; 0.33,0; 0.4,0.045; 0.8,0.54; ...
    1,0.9; 1.01,0.25; 2,0; 2.4 0];
c = interp1(FSTTB_tab(:,1),FSTTB_tab(:,2),tau);
end
function [c] = FSOTB(tau)
FSOTB_tab = [ 0, 0; 0.33, 0; 0.4, 0; 0.8, 0; 1, 0; 1.01,0.75; 2,1; 2.4,1];
c = interp1(FSOTB_tab(:,1),FSOTB_tab(:,2),tau);
end
%% fertilizer functions % added 240,0.7 because if DAY>200 matlab gives NAN when interpolating since its not in range
function [c] = FERTAB(tau)
FERTAB_tab = [0,0;99,0;100,10;101,0;124,0;125,5;126,0;200,0 ;240,0];
c = interp1(FERTAB_tab(:,1),FERTAB_tab(:,2),tau);
end
function [c] = NRFTAB(tau)
NRFTAB_tab = [0,0.7; 100,0.7; 125,0.7; 150,0.7;200,0.7; 240,0];
c = interp1(NRFTAB_tab(:,1),NRFTAB_tab(:,2),tau);
end
%% leaf death
function [RDR] = DLA(tau,dtau,LAI)
global LAI_cr DRSH0 DRDV0 TSUMAN

% death due to ageing
if tau >= TSUMAN
    DRDV = interp1(DRDV0(:,1),DRDV0(:,2),dtau);
else
    DRDV = 0;
end

% death due to shading
DRSH = max(0, min(DRSH0, DRSH0*(LAI-LAI_cr)/LAI_cr));

RDR = max(DRSH,DRDV);
end
function [DAYL] = ASTRO(DOY,LAT)
SINLAT = sin(pi*LAT/180);
COSLAT = cos(pi*LAT/180);
SINDCM = sin(pi*23.45/100);
SINDEC = -SINDCM * cos(2*pi*(DOY+10)/365);
COSDEC = sqrt(1-SINDEC*SINDEC);
A = SINLAT*SINDEC;
B = COSLAT*COSDEC;
DAYL = 12*(1+(2/pi)*asin(A/B));
return;
end
function [DVS,DVS1,DVS2] = SUBDVS(TIME,DOYEM,TSUM,TSUMAN,TSUMMT,DVS1,DVS2)
if(TIME < DOYEM)
    DVS1 = TSUM/TSUMAN;
elseif(TIME >= DOYEM && TSUM <= TSUMAN)
    DVS1 = TSUM/TSUMAN;
else
    DVS2 = (TSUM-TSUMAN)/TSUMMT;
end
DVS = DVS1+DVS2;
return
end

%% LAI growth
function [glai] = GLA(time,doyem,dteff,laii,rgrl,delt,sla,lai,glv,nlai,wc,wcwp,dvs,tranrf,nni)
% ---------------------------------------------------------------------*
%  SUBROUTINE GLA                                                      *
%  Purpose: This subroutine computes daily increase of leaf area index *
%           (ha leaf/ ha ground/ d)                                    *
% ---------------------------------------------------------------------*
%---- Growth during maturation stage:
glai = sla .* glv;

%---- Growth during juvenile stage:
if((dvs<0.2)&&(lai<0.75))
    glai = (lai .*(exp(rgrl .* dteff .* delt) - 1.) ./ delt) .* tranrf*exp(-nlai*(1-nni));
end

%---- Growth at day of seedling emergence:
if((time>=doyem)&&(lai==0.)&&(wc>wcwp))
    glai = laii ./ delt;
end

%---- Growth before seedling emergence:
if(time < doyem)
    glai = 0.;
end

return;
end
function [pevap,ptran]=PENMAN(davtmp,vp,dtr,lai,wn)
% ---------------------------------------------------------------------*
%  SUBROUTINE PENMAN                                                   *
%  Purpose: Computation of the PENMAN EQUATION                         *
% ---------------------------------------------------------------------*
dtrjm2 = dtr .* 1.0e6;
boltzm = 5.668e-8;
lhvap  = 2.4e6;
psych  = 0.067;

bbrad  = boltzm .*(davtmp+273.).^4 .* 86400.;
svp    = 0.611 .* exp(17.4 .* davtmp ./(davtmp + 239.));
slope  = 4158.6 .* svp ./(davtmp + 239.).^2;
rlwn   = bbrad .* max(0.,0.55.*(1.-vp./svp));
nrads  = dtrjm2 .*(1.-0.15) - rlwn;
nradc  = dtrjm2 .*(1.-0.25) - rlwn;
penmrs = nrads .* slope./(slope+psych);
penmrc = nradc .* slope./(slope+psych);

wdf    = 2.63 .*(1.0 + 0.54 .* wn);
penmd  = lhvap .* wdf .*(svp-vp) .* psych./(slope+psych);

pevap  =     exp(-0.5.*lai)  .*(penmrs + penmd) ./ lhvap;
pevap = max(0,pevap);
ptran  =(1.-exp(-0.5.*lai)) .*(penmrc + penmd) ./ lhvap;
ptran  = max( 0., ptran );

return;
end
function [EVAP,TRAN,WCCR,DSLR]=EVAPTR(PEVAP,PTRAN,rootd,WA,WCAD,WCWP,WCFC,WCWET,WCST,TRANCO,DELT,WMFAC,RAIN,DSLR)
% ---------------------------------------------------------------------*
%  SUBROUTINE EVAPTR                                                   *
%  Purpose: To compute actual rates of evaporation and transpiration   *
% ---------------------------------------------------------------------*
WC   = 0.001 .* WA   ./ rootd;
WAAD = 1000. .* WCAD .* rootd;
WAFC = 1000. .* WCFC .* rootd;

if(RAIN>= 0.5)
    EVS = PEVAP;
    DSLR = 1;
else
    DSLR = DSLR+1;
    EVSMXT = PEVAP*(sqrt(DSLR)-sqrt(DSLR-1));
    EVS = min(PEVAP,EVSMXT+RAIN);
end

WCCR = WCWP+max(0.01,PTRAN/(PTRAN+TRANCO)*(WCFC-WCWP));

if(WMFAC>= 1)
    if(WC>WCCR)
        FR = 1;
    else
        if((WC-WCWP)/(WCCR-WCWP)<=1 && (WC-WCWP)/(WCCR-WCWP)>= 0)
            FR = (WC-WCWP)/(WCCR-WCWP);
        elseif((WC-WCWP)/(WCCR-WCWP)<0)
            FR = 0;
        else
            FR = 1;
        end
    end
else
    if(WC>WCCR)
        if((WCST-WC)/(WCST-WCWET)<=1 && (WCST-WC)/(WCST-WCWET)>= 0)
            FR = (WCST-WC)/(WCST-WCWET);
        elseif((WCST-WC)/(WCST-WCWET)<0)
            FR = 0;
        else
            FR = 1;
        end
    else
        if((WC-WCST)/(WCCR-WCWP)<=1 && (WC-WCST)/(WCCR-WCWP)>= 0)
            FR = (WC-WCST)/(WCCR-WCWP);
        elseif((WC-WCST)/(WCCR-WCWP)<0)
            FR = 0;
        else
            FR = 1;
        end
        
    end
    
end
TRAN = PTRAN*FR;
AVAILF = min(1,((WA-WAAD)/DELT)/(EVS+TRAN));
EVAP = EVS*AVAILF;
TRAN = TRAN*AVAILF;

return;
end

%% Drainage, runoff and irrigation
function [DRAIN,RUNOFF,IRRIG]=DRUNIR(RAIN,EVAP,TRAN,IRRIGF,DRATE,DELT,WA,ROOTD,WCFC,WCST,WMFAC)
WC   = 0.001 .* WA   ./ ROOTD;
WAFC = 1000. .* WCFC .* ROOTD;
WAST = 1000. .* WCST .* ROOTD;
if((WA-WAFC)/DELT +  (RAIN - EVAP - TRAN)>= 0 && (WA-WAFC)/DELT +  (RAIN - EVAP - TRAN)<=DRATE)
    DRAIN = (WA-WAFC)/DELT +  (RAIN - EVAP - TRAN);
elseif((WA-WAFC)/DELT +  (RAIN - EVAP - TRAN)<0)
    DRAIN = 0;
else
    DRAIN = DRATE;
end
RUNOFF = max( 0., (WA-WAST)/DELT +  (RAIN - EVAP - TRAN - DRAIN));
if(WMFAC >= 1)    
    IRRIG = IRRIGF*max(0,(WAST-WA)/DELT - (RAIN - EVAP- TRAN - DRAIN - RUNOFF));
else
    IRRIG = IRRIGF*max(0,(WAFC-WA)/DELT - (RAIN - EVAP- TRAN - DRAIN - RUNOFF));
end
return;
end
function [FRT,FLV,FST,FSO] = SUBPAR(NPART,TRANRF,NNI,FRTWET,FLVT,FSTT,FSOT)
if(TRANRF < NNI)
    FRTMOD = max( 1, 1/(TRANRF+0.5));
    FRT    = FRTWET * FRTMOD;
    FSHMOD = (1.-FRT) / (1.-FRT/FRTMOD);
    FLV    = FLVT * FSHMOD;
    FST    = FSTT * FSHMOD;
    FSO    = FSOT * FSHMOD;
else
    FLVMOD = exp(-NPART* (1.0-NNI));
    FLV    = FLVT * FLVMOD;
    MODIF  = (1.-FLV)/(1.-(FLV/FLVMOD));
    FST    = FSTT *  MODIF;
    FRT    = FRTWET* MODIF;
    FSO    = FSOT *  MODIF;
end
return;
end
function [PARINT,GTOTAL] = GROWTH(TIME,DOYEM,DTR,K,NLUE,LAI,LUE,TRANRF,NNI)
if((TIME-DOYEM)>=0)
    PARINT = 0.5 * DTR * (1.- exp(-K*LAI));
else
    PARINT = 0;
end
if(TRANRF <= NNI)
    GTOTAL = LUE * PARINT * TRANRF;
else
    GTOTAL = LUE * PARINT *exp(-NLUE*(1-NNI));
end
return;
end
function [RWLVG,RWRT,RWST,RWSO] = RELGR(TIME,DOYEM,EMERG,WLVGI,WRTLI,WSTI,WSOI,GTOTAL,FLV,FRT,FST,FSO,DLV,DRRT,DELT)
if(TIME>=DOYEM && EMERG == 1)
    RWLVG = GTOTAL * FLV - DLV;
    RWRT  = GTOTAL * FRT - DRRT;
    RWST  = GTOTAL * FST;
    RWSO  = GTOTAL * FSO;
else
    RWLVG = 0;
    RWRT  = 0;
    RWST  = 0;
    RWSO  = 0;
end
return;
end
function [RNULV,RNUST,RNURT] = RNUSUB(TIME,DOYEM,EMERG,NDEML,NDEMS,NDEMR,NUPTR,NDEMTO,ANLVI,ANSTI,ANRTI,DELT)
if(TIME>=DOYEM && EMERG == 1)
    if(NDEMTO == 0)
        RNULV = NDEML* NUPTR;
        RNUST = NDEMS* NUPTR;
        RNURT = NDEMR* NUPTR;
    else
        RNULV = (NDEML / (NDEMTO))* NUPTR;
        RNUST = (NDEMS / (NDEMTO))* NUPTR;
        RNURT = (NDEMR / (NDEMTO))* NUPTR;
    end
else
    RNULV = 0;
    RNUST = 0;
    RNURT = 0;
end

return;
end
function [NOPTLV,NOPTST] = NOPTM(FRNX,NMAXLV,NMAXST)
NOPTLV= FRNX * NMAXLV;
NOPTST= FRNX * NMAXST;
return;
end
function [NDEML,NDEMS,NDEMR,NDEMSO] = NDEMND(NMAXLV,NMAXST,NMAXRT,NMAXSO,WLVG,WST,WRT,WSO,ANLV,ANST,ANRT,ANSO,TCNT,DELT)
NDEML  =  max(NMAXLV*WLVG -  ANLV, 0.);
NDEMS  =  max(NMAXST*WST  - ANST, 0.);
NDEMR  =  max(NMAXRT*WRT  - ANRT, 0.);
NDEMSO =  max(NMAXSO*WSO  - ANSO, 0.)/TCNT;
return;
end
function [ATNLV,ATNST,ATNRT,ATN] = NTRLOC(ANLV,ANST,ANRT,WLVG,WST,WRT,RNFLV,RNFST,RNFRT,FNTRT)
ATNLV = max (0. , ANLV-WLVG*RNFLV);
ATNST = max (0. , ANST-WST*RNFST);
ATNRT = min((ATNLV + ATNST) * FNTRT, ANRT-WRT*RNFRT);
ATN   = ATNLV +  ATNST + ATNRT;
return;
end
function [RNTLV,RNTST,RNTRT] = NTRANS(RNSO,ATNLV,ATNST,ATNRT,ATN)
if(ATN == 0)
    RNTLV = RNSO*ATNLV;
    RNTST = RNSO*ATNST;
    RNTRT = RNSO*ATNRT;
else
    RNTLV = RNSO*ATNLV/ATN;
    RNTST = RNSO*ATNST/ATN;
    RNTRT = RNSO*ATNRT/ATN;
end
return;
end
function [DRRT,RNLDLV,RNLDRT] = RNLD(DVS,WRT,RDRRT,RNFLV,DLV,RNFRT,DVSDR)
if(DVS < DVSDR)
    DRRT =0;
else
    DRRT = WRT*RDRRT;
end
RNLDLV = RNFLV*DLV;
RNLDRT = RNFRT*DRRT;
return;
end
function [NNI] = NNINDX(TIME,DOYEM,EMERG,NFGMR,NRMR,NOPTMR)
TINY = 0.001;
if(TIME >= DOYEM && EMERG == 1)
    if( (NFGMR-NRMR)/(NOPTMR-NRMR)>=TINY&& (NFGMR-NRMR)/(NOPTMR-NRMR)<=1.0 )
        NNI = (NFGMR-NRMR)/(NOPTMR-NRMR);
    elseif((NFGMR-NRMR)/(NOPTMR-NRMR)<TINY)
        NNI = TINY;
    else
        NNI = 1.0;
    end
else
    NNI = 0;
end

return;
end
function [RDRDV,RDRSH,RDR,DLV,DLVS,DLVNS,DLAIS,DLAINS,DLAI] = DEATHL(TIME,DOYEM,TSUM,TSUMAG,RDRTMP,RDRSHM,LAI,LAICR,WLVG,RDRNS,NNI,SLA)
if(TSUM < TSUMAG)
    RDRDV =0;
else
    RDRDV = RDRTMP;
end
RDRSH = max(0.,RDRSHM*(LAI-LAICR)/LAICR);
RDR = max(RDRDV,RDRSH);

if(TIME >= DOYEM && NNI < 1)
    DLVNS = WLVG*RDRNS*(1-NNI);
    DLAINS = DLVNS*SLA;
else
    DLVNS = 0;
    DLAINS = 0;
end
DLVS = WLVG*RDR;
DLAIS = LAI*RDR;
DLV = DLVS + DLVNS;
DLAI = DLAIS + DLAINS;
return;
end













