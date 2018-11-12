global SLAC TBASE NMAXSO DELT WMFAC RDRNS DVSNT DVSNLT WCSUBS DVSDR LSNR FRNX LRNR WCI NFRLVI NFRSTI NFRRTI RGRL LAICR TSUMAN TSUMMT TSUMAG RDRSHM LUE K ROOTDM RRDMAX WCAD WCWET WCST WCWP WCFC TRANCO DRATE IRRIGF TCNT RDRRT RNFLV RNFST RNFRT FNTRT NLUE NLAI NSLA NPART
% specific leaf area constant
SLAC = 0.022;
% base temperature for spring wheat crop
TBASE = 0;
% relative growth rate of LAI at the xponential growth phase (oC d)^-1,
% critical LAI above which mutual shading of leaves occur, and the maximum
% relative death rate of leaves due to shading.
RGRL = 0.009; LAICR = 4; RDRSHM = 0.03;
% light use efficiency and extinction coefficient
LUE = 2.8; K = 0.6;
% maximum root depth and maximum rate of increase in rooting depth (m d-1)
% for a rice crop
ROOTDM = 1.2; RRDMAX = 0.012;
% soil hydraulic properties
WCAD = 0.1; WCWP = 0.2; WCFC = 0.4; WCWET = 0.45; WCST = 0.5;
% transpiration constant (mm/day) indicating the level of drought tolerance
% of the wheat crop
TRANCO = 8;
% maximum drainage rate of the soil (mm/day) and irrigation factor (1 = yes, 0 = no)
DRATE = 30; IRRIGF = 1;
% time constant (days) for N translocation
TCNT = 10;
% relative death rate of roots
RDRRT = 0.03;
% residual N concentration in leaves, stem and roots.
RNFLV = 0.004; RNFST = 0.002; RNFRT = 0.002;
% nitrogen translocated from roots as a fraction of the total amount of
% nitrogen translocated from leaves and stem.
FNTRT = 0.15;
% extinction coefficient for nitrogen distribution down the canopy
NLUE = 0.2;
% coefficient for the effect of N stress on LAI reduction (during juvenile stage)
NLAI = 1;
% coefficient for the effect of N stress on SLA reduction
NSLA = 1;
% coefficient for the effect of N stress on leaf biomass reduction
NPART = 1;
% initial water content in (water/(cm3 soil)
WCI = 0.4; % old 0.4
% temperature sum for anthesis, and maturity and ageing of leaves
TSUMAN = 800; TSUMMT = 1030; TSUMAG = 800;
% initial fraction of N (g N g-1 DM) in leaves, stem and roots.
NFRLVI = 0.06;  NFRSTI = 0.03; NFRRTI = 0.03;
% used for N concentration in th leaves, from which stem and roots are
% derived as function of DVS
NMAXSO = 0.0165;LSNR = 0.50; LRNR = 0.50;
% optimal N concentration as the fraction of maximum N concentration
FRNX = 0.5;
% Nitrogen limiting factor at low moisture conditions
DVSNLT = 1.0;
% development stage above which death of leaves and roots starts
DVSDR = 1.0;
% as irrigated up to field capacity (WMFAC = 0), as irrigated up to
% saturation (WMFAC = 1)
WMFAC = 0.0;
% constant of exploartion of water in soil when roots grow downwards
WCSUBS= 0.30;
% increase time simulation
DELT = 1;
% relative death of leaves due to N stress
RDRNS = 0.03;
% factor N supply to storage organs
DVSNT = 0.8;
% dummy variable (nested in FORTRAN)
if(exist('m','var') == 0)
    DSLR = 1; DVS1 = 0; DVS2 = 0; 
else
    DSLR = ones(m,1);
    DVS1 = zeros(m,1);
    DVS2 = zeros(m,1);
    DSLR_r = DSLR;
    DVS1_r = DVS1;
    DVS2_r = DVS2;
end

