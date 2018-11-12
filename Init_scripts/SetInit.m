% initial development stage
TSUMI = 0;
DVSI = 0;
% Initial root depth (m)
ROOTDI = 0.1;
% initial amount of water present in the rooted depth at the start of the
% calculations, based on the initial water content (mm)
WAI = 1000*ROOTDI*WCI;
% initial weigth of leaves, roots stem, and storage organs and
% transplanting (g/m2)
WLVGI = 2.4; WSTI  = 0.0; WRTLI = 3.6; WSOI  = 0.0;
% initial LAI
SLACFI = 1; % SLACF(DVSI);
ISLA = SLAC * SLACFI;
LAII = WLVGI * ISLA;
% Initial amount of N (g/cm2) in leaves, stem, roots and storage organs.
ANLVI = NFRLVI* WLVGI;
ANSTI = NFRSTI* WSTI;
ANRTI = NFRRTI* WRTLI;
ANSOI = 0.;
% allocate initial conditions to state variables
if(exist('m','var') == 0)
    TSUM(:,1) = TSUMI;
    WLVG(:,1) = WLVGI;
    DVS(:,1) = DVSI;
    ROOTD(:,1) = ROOTDI;
    WA(:,1) = WAI;
    WC(:,1) = WCI;
    WCCR(:,1) = WCWP+0.01;
    LAI(:,1) = LAII;
    WLVG(:,1) = WLVGI;
    WRT(:,1) = WRTLI;
    WST(:,1) = WSTI;
    WSO(:,1) = WSOI;
    ANLV(:,1) = ANLVI;
    ANST(:,1) = ANSTI;
    ANRT(:,1) = ANRTI;
    ANSO(:,1) = ANSOI;
else
    TSUM_r(:,1) = TSUMI;
    WLVG_r(:,1) = WLVGI;
    DVS_r(:,1) = DVSI;
    ROOTD_r(:,1) = ROOTDI;
    WA_r(:,1) = WAI;
    WC_r(:,1) = WCI;
    WCCR_r(:,1) = WCWP+0.01;
    LAI_r(:,1) = LAII;
    WLVG_r(:,1) = WLVGI;
    WRT_r(:,1) = WRTLI;
    WST_r(:,1) = WSTI;
    WSO_r(:,1) = WSOI;
    ANLV_r(:,1) = ANLVI;
    ANST_r(:,1) = ANSTI;
    ANRT_r(:,1) = ANRTI;
    ANSO_r(:,1) = ANSOI;
end