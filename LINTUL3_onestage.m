function [TSUM_n,ROOTD_n,WA_n,WC_n,WCCR_n,TEXPLO_n,TEVAP_n,TTRAN_n,TRUNOF_n,TIRRIG_n,TDRAIN_n,TRAIN_n,DVS_n,NNI_n,SLA_n,LAI_n,NDEMTO_n,TNSOIL_n,NUPTT_n,ANLV_n,ANST_n,ANRT_n,ANSO_n,NLOSSL_n,NLOSSR_n,WLVG_n,WLVD_n,WST_n,WSO_n,WRT_n,TAGBM_n,NTAC_n,LUECAL_n,CUMPAR_n, GTSUM_n,WDRT_n,DSLR_n,DVS1_n,DVS2_n]...
    = LINTUL3_onestage(DOY,LAI,TSUM,DTEFF,WC,ROOTD,RAIN,DAVTMP,VP,DTR,WN,TEXPLO,WA,TEVAP,TTRAN,TRUNOF,TIRRIG,TDRAIN,TRAIN,WLVG,WST,ANLV,ANST,ANRT,ANSO,WRT,WSO,TNSOIL,NUPTT,NLOSSL,NLOSSR,WLVD,WDRT,GTSUM,CUMPAR,DSLR,DVS1,DVS2,man_irrig,irrig_fact,man_fert,fert_fact)
    global DOYEM WCWP TSUMAN RRDMAX ROOTDM WCSUBS LSNR LRNR RNFLV RNFST SLAC NSLA DVSNLT DELT DVSNT LAII TCNT 
    %% --------- emergence, temperature sum --------- 
    EMERG = ((DOY>=(DOYEM)&& WC>WCWP) && -LAI<0 );
    % photosyntatic effect, latitude amersvoort 571.460
    DAYL = ASTRO(DOY,571.460); % weatherstation 6 netherlands
    if((TSUM-TSUMAN)<0)
        PHOTPF = PHOTTB(DAYL);
    else
        PHOTPF = 1;
    end
    RTSUM = DTEFF*EMERG; 
    RTSUMP = RTSUM * PHOTPF;
    TSUM_n = TSUM + RTSUMP; 
    %% --------- Rooth Depth ---------
    RROOTD = min(RRDMAX*(WC >= WCWP)*EMERG,ROOTDM-ROOTD);
    ROOTD_n = ROOTD + RROOTD;
    %% --------- Soil moisture balance --------- 
    [PEVAP,PTRAN] = PENMAN(DAVTMP,VP,DTR,LAI,WN);    
    [EVAP, TRAN, WCCR_n,DSLR_n] = EVAPTR(PEVAP,PTRAN,ROOTD,WA,RAIN,DSLR); % DSLR global maybe   
    [DRAIN,RUNOFF,IRRIG] = DRUNIR(RAIN,EVAP,TRAN,WA,ROOTD);
    % TRANRF verified
    if(PTRAN <= eps)
        TRANRF = TRAN;
    else
        TRANRF = TRAN/PTRAN; 
    end    
    % update water levels
    EXPLOR = 1000*RROOTD*WCSUBS;
    if(irrig_fact == 1)
        RWA = (RAIN+EXPLOR+IRRIG) - (RUNOFF+TRAN+EVAP+DRAIN);  % optimal rainfed
    else
        RWA = (RAIN+EXPLOR+man_irrig) - (RUNOFF+TRAN+EVAP+DRAIN);  % manual rainfed
    end
    WA_n = WA+RWA;
    % water content in rootzone
    WC_n = WA_n*0.001/ROOTD_n;
    % total exploration
    TEXPLO_n = TEXPLO + EXPLOR;
    % crop transipartion and soil evaporation
    TEVAP_n = TEVAP + EVAP;
    TTRAN_n = TTRAN + TRAN;
    % other water balance integrals
    TRUNOF_n = TRUNOF + RUNOFF;
    TIRRIG_n = TIRRIG + IRRIG;
    TDRAIN_n = TDRAIN + DRAIN;
    TRAIN_n = TRAIN + RAIN;
    % update water balance
    % WATBAL(k) = WA(k)-WAI-TRAIN(k)-TEXPLO(k)-TIRRIG(k)+TRUNOF(k)+TTRAN(k)+TEVAP(k)+TDRAIN(k);
    %% --------- Calculation of Develement stage  DVS and update NNI ---------
    [DVS_n,DVS1_n,DVS2_n] = SUBDVS(DOY,TSUM,DVS1,DVS2); % DVS1/DVS2 global?
    % maximum N concentration in the leaves, from which the values of the
    % stem and roots are derived, as a function of development stage.
    NMAXLV = NMXLV(DVS_n);
    NMAXST = LSNR * NMAXLV;
    NMAXRT = LRNR * NMAXLV;
    % maximum nitrogen concentration of stem and leaves
    [NOPTLV,NOPTST] = NOPTM(NMAXLV,NMAXST);
    % total vegatative biomass
    TBGMR = WLVG+WST;
    % maximum N content in the plant
    NOPTS = NOPTST* WST;
    NOPTL = NOPTLV* WLVG;  
    if(TBGMR == 0)
        NOPTMR = (NOPTL+NOPTS);
    else
        NOPTMR = (NOPTL+NOPTS)/TBGMR;
    end
    NUPGMR = ANLV+ANST;
    if(TBGMR == 0)
        NRMR = (WLVG*RNFLV+WST*RNFST);
        NFGMR  = NUPGMR;
    else
        NRMR = (WLVG*RNFLV+WST*RNFST)/(TBGMR);
        NFGMR = NUPGMR/TBGMR;
    end
    NNI_n = NNINDX(DOY,EMERG,NFGMR,NRMR,NOPTMR);
    SLA_n = SLAC * SLACF(DVS_n)*exp(-NSLA*(1-NNI_n));
    
    %% --------- soil nitrogen supply ---------
    RDRTMP = RDRT(DAVTMP);
    % growth rate and dry matter production of plant organs
    FRTWET = FRTTB(DVS_n); 
    FLVT = FLVTB(DVS_n); 
    FSTT = FSTTB(DVS_n); 
    FSOT = FSOTB(DVS_n); 
    [FRT,FLV,FST,FSO] = SUBPAR(TRANRF,NNI_n,FRTWET,FLVT,FSTT,FSOT);
    %FRACT = FLV+FST+FSO+FRT; 
    
    %% calculate leaf growth
    [~,GTOTAL] = GROWTH(DOY,DTR,LAI,TRANRF,NNI_n);
    GLV = FLV * GTOTAL;      
    GLAI = GLA(DOY,DTEFF,LAII,SLA_n,LAI,GLV,WC,DVS_n,TRANRF,NNI_n); % make LAII global
    [~,~,~,DLV,~,~,~,~,DLAI] = DEATHL(DOY,TSUM,RDRTMP,LAI,WLVG,NNI_n,SLA_n);
    RLAI = GLAI - DLAI;        
    LAI_n = LAI + RLAI; 
    
    %% --------- soil nitrogen supply continueed ---------
    [DRRT,RNLDLV,RNLDRT] = RNLD(DVS_n,WRT,DLV);
    [RWLVG,RWRT,RWST,RWSO] = RELGR(DOY,EMERG,GTOTAL,FLV,FRT,FST,FSO,DLV,DRRT);
    [NDEML,NDEMS,NDEMR,NDEMSO] = NDEMND(NMAXLV,NMAXST,NMAXRT,WLVG,WST,WRT,WSO,ANLV,ANST,ANRT,ANSO);
    NDEMTO_n = max(0,NDEML+NDEMS+NDEMR);   
    % Nitrogren limiting factor at low moisture conditions
    if(DVS_n<DVSNLT && WC>= WCWP)
        NLIMIT = 1;
    else
        NLIMIT = 0;
    end   
    % Soil N supply
    RTMIN  = 0.10 * EMERG * NLIMIT;    
    % fertilizer application use LINTUL 3 example or manual fertilizer
    if(fert_fact == 1)
        NRF = NRFTAB(DOY);
        FERTN = FERTAB(DOY);
        FERTNS = FERTN * NRF; 
    else 
        FERTNS = man_fert;
    end
    % total nutrient uptake
    if(DOY<DOYEM)
        NUPTR = 0;
    else
        NUPTR = max(0., min (NDEMTO_n, TNSOIL))* NLIMIT/ DELT;
    end 
    [RNULV,RNUST,RNURT] = RNUSUB(DOY,EMERG,NDEML,NDEMS,NDEMR,NUPTR,NDEMTO_n);
    [ATNLV,ATNST,ATNRT,ATN] = NTRLOC(ANLV,ANST,ANRT,WLVG,WST,WRT);
    if(DVS_n<DVSNT)
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
    RNSOIL = FERTNS/DELT -NUPTR + RTMIN;
    % Amount of inorganic N in soil as function of fertizlier
    TNSOIL_n = TNSOIL + RNSOIL;
    % Total uptake of N over time (g N m-2)
    NUPTT_n = NUPTT + NUPTR; 
    % actual N concent of organs
    ANLV_n = ANLV + RNLV;
    ANST_n = ANST + RNST;
    ANRT_n = ANRT + RNRT;
    ANSO_n = ANSO + RNSO;
    % Nitrogen Balance
    NLOSSL_n = RNLDLV + NLOSSL;
    NLOSSR_n = NLOSSR + RNLDRT;
    % NBALAN = NUPTT+(ANLVI+ANSTI+ANRTI+ANSOI)-(ANLV+ANST+ANRT+ANSO+NLOSSL+NLOSSR); 
    % weight of green leaves, dead leaves, stem, storage orangs and roots
    WLVG_n = WLVG+RWLVG;
    WLVD_n = WLVD+DLV;
    WST_n = WST+RWST; 
    WSO_n = WSO+RWSO; 
    WRT_n = WRT+RWRT;
    % total leaf weigth
    WLV = WLVG+WLVD;
    % total above ground biomass
    TAGBM_n = WLV+WST+WSO; 
    % Water Nitrogren stress factor
    RNW = min(TRANRF,NNI_n);
    % Biomass carbon balance
    GTSUM_n = GTSUM+GTOTAL;
    WDRT_n = WDRT + DRRT;
    % CBALAN = GTSUM+(WRTLI+WLVGI+WSTI+WSOI)-(WLV+WST+WSO+WRT+WDRT);    % include initials in globals 
    % For calculation of LUE
    PAR = 0.5*DTR;
    LUECAL_n = GTOTAL/PAR;
    CUMPAR_n = CUMPAR + PAR;
    NTAG = ANLV+ANST+ANSO;
    if(NTAG == 0)
        NTAC_n = NTAG;
    else
        NTAC_n = NTAG/TAGBM_n;
    end
end

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

%% globals included subroutines FORTRAN
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
function [DVS,DVS1,DVS2] = SUBDVS(TIME,TSUM,DVS1,DVS2)
global TSUMAN TSUMMT DOYEM 

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
function [glai] = GLA(time,dteff,laii,sla,lai,glv,wc,dvs,tranrf,nni)
% ---------------------------------------------------------------------*
%  SUBROUTINE GLA                                                      *
%  Purpose: This subroutine computes daily increase of leaf area index *
%           (ha leaf/ ha ground/ d)                                    *
% ---------------------------------------------------------------------*
%---- Growth during maturation stage:
global RGRL DELT DOYEM NLAI WCWP
glai = sla .* glv;

%---- Growth during juvenile stage:
if((dvs<0.2)&&(lai<0.75))
    glai = (lai .*(exp(RGRL .* dteff .* DELT) - 1.) ./ DELT) .* tranrf*exp(-NLAI*(1-nni));
end

%---- Growth at day of seedling emergence:
if((time>= DOYEM)&&(lai==0.)&&(wc>WCWP))
    glai = laii ./ DELT;
end

%---- Growth before seedling emergence:
if(time < DOYEM)
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
function [EVAP,TRAN,WCCR,DSLR]=EVAPTR(PEVAP,PTRAN,rootd,WA,RAIN,DSLR)
% ---------------------------------------------------------------------*
%  SUBROUTINE EVAPTR                                                   *
%  Purpose: To compute actual rates of evaporation and transpiration   *
% ---------------------------------------------------------------------*
global WCAD WCWP WCFC WCWET WCST TRANCO DELT WMFAC 
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
function [DRAIN,RUNOFF,IRRIG]=DRUNIR(RAIN,EVAP,TRAN,WA,ROOTD)
global IRRIGF DRATE DELT WCFC WCST WMFAC
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
function [FRT,FLV,FST,FSO] = SUBPAR(TRANRF,NNI,FRTWET,FLVT,FSTT,FSOT)
global NPART
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
function [PARINT,GTOTAL] = GROWTH(TIME,DTR,LAI,TRANRF,NNI)
global K NLUE DOYEM LUE
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
function [RWLVG,RWRT,RWST,RWSO] = RELGR(TIME,EMERG,GTOTAL,FLV,FRT,FST,FSO,DLV,DRRT)
global DOYEM
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
function [RNULV,RNUST,RNURT] = RNUSUB(TIME,EMERG,NDEML,NDEMS,NDEMR,NUPTR,NDEMTO)
global DOYEM
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
function [NOPTLV,NOPTST] = NOPTM(NMAXLV,NMAXST)
global FRNX
NOPTLV= FRNX * NMAXLV;
NOPTST= FRNX * NMAXST;
return;
end
function [NDEML,NDEMS,NDEMR,NDEMSO] = NDEMND(NMAXLV,NMAXST,NMAXRT,WLVG,WST,WRT,WSO,ANLV,ANST,ANRT,ANSO)
global NMAXSO TCNT
NDEML  =  max(NMAXLV*WLVG -  ANLV, 0.);
NDEMS  =  max(NMAXST*WST  - ANST, 0.);
NDEMR  =  max(NMAXRT*WRT  - ANRT, 0.);
NDEMSO =  max(NMAXSO*WSO  - ANSO, 0.)/TCNT;
return;
end
function [ATNLV,ATNST,ATNRT,ATN] = NTRLOC(ANLV,ANST,ANRT,WLVG,WST,WRT)
global RNFLV RNFST RNFRT FNTRT
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
function [DRRT,RNLDLV,RNLDRT] = RNLD(DVS,WRT,DLV)
global RDRRT RNFLV DVSDR RNFRT
if(DVS < DVSDR)
    DRRT =0;
else
    DRRT = WRT*RDRRT;
end
RNLDLV = RNFLV*DLV;
RNLDRT = RNFRT*DRRT;
return;
end
function [NNI] = NNINDX(TIME,EMERG,NFGMR,NRMR,NOPTMR)
global DOYEM
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
function [RDRDV,RDRSH,RDR,DLV,DLVS,DLVNS,DLAIS,DLAINS,DLAI] = DEATHL(TIME,TSUM,RDRTMP,LAI,WLVG,NNI,SLA)
global DOYEM TSUMAG RDRSHM LAICR RDRNS
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


