%% Lintul 3 Version
addpath('Init_scripts')
clear variables; close all; clc;
%% General structure of DHMARA
m = 6; % number of clients
p = 1; % number of products
q = 3; % number of subtime quantizations
Cmax = [30 30 0;0 0 10]; % agent max capacity matrix four agents that can deliver water and 1 agents that can deliver fertilizer
[r,n] = size(Cmax); % number of resources (r) and number of agents (n)

% Weighting vectors
pi    = 100;%rand(p,1); % price per product vector
rho   = [10; 100];%rand(r,1); % cost per resource vector (change second term)
sigma = 10*ones(n,1); % cost per operation vector

U = [0 10/3 20/3 10 0 0 0; [zeros(1,4) linspace(4,10,3)]]; % r*s
s = size(U,2); % number of bids per client
Ucheck = reshape(U,[],1)*ones(1,m); %(r*s)*m
decision_variables = m*(s+q*(n+r));

%% Optimization constraints
delt = binvar(s,m,'full');
Atau = binvar(n,m,q+1,'full');
y = binvar(n,1,q+1,'full');
Atau_plus = [Atau y];
Dtau = sdpvar(r,m,q,'full');
H_invent = sdpvar(r,n,q+1,'full');
C_in = sdpvar(r,n,q,'full');
C_bid_selection = [delt'*ones(s,1) == 1];
C_1fieldpagent = [];
C_1agentpfield = [];
C_no_free_deli = [];
C_agent_cap    = [];
C_initconstraints = [];
H_inventconstraints = [];
E = GetE_gui_grid(m,3);
F = GetF_gui_grid(m,3);
E = [E F; F' 0];
A_init = [zeros(n,m) ones(n,1)];
C_init = Cmax;
for t = 1:q
    C_1fieldpagent = [C_1fieldpagent Atau_plus(:,:,t)*ones(m+1,1) == 1  ];
    C_1agentpfield = [C_1agentpfield Atau_plus(:,1:m,t)'*ones(n,1) <= 1 Atau_plus(:,m+1,t)'*ones(n,1) <= n];
    C_no_free_deli = [C_no_free_deli Dtau(:,:,t)*(eye(m)-diag(Atau_plus(:,1:m,t)'*ones(n,1) )) == 0];
    C_agent_cap = [C_agent_cap Dtau(:,:,t)*Atau_plus(:,1:m,t)' <= Cmax];
    if(t==1)
        %C_initconstraints = [C_in(:,:,t) == Cmax-H_invent(:,:,q+1)];
        %H_inventconstraints = [H_invent(:,:,t) == H_invent(:,:,q+1)-Dtau(:,:,t)*Atau_plus(:,1:m,t)'+C_in(:,:,t)*diag((ones(n,1)-y(:,:,q+1)).*y(:,:,t))>= 0];
        H_invent(:,:,t) = H_invent(:,:,q+1)-Dtau(:,:,t)*Atau_plus(:,1:m,t)'+(Cmax-H_invent(:,:,q+1))*diag((ones(n,1)-y(:,:,q+1)).*y(:,:,t));
        H_inventconstraints = [H_invent(:,:,t) >= 0];
    else
        %C_initconstraints = [C_initconstraints C_in(:,:,t) == Cmax-H_invent(:,:,t-1)];
        %H_inventconstraints = [H_inventconstraints H_invent(:,:,t) == H_invent(:,:,t-1)-Dtau(:,:,t)*Atau_plus(:,1:m,t)'+C_in(:,:,t)*diag((ones(n,1)-y(:,:,t-1)).*y(:,:,t)) >= 0];
        H_invent(:,:,t) = H_invent(:,:,t-1)-Dtau(:,:,t)*Atau_plus(:,1:m,t)'+(Cmax-H_invent(:,:,t-1))*diag((ones(n,1)-y(:,:,t-1)).*y(:,:,t));
        H_inventconstraints = [H_inventconstraints H_invent(:,:,t)>= 0];
    end
    
end
C_u_lo  = [];
C_u_hi  = [];
DUcheck = kron(delt,ones(r,1)).*Ucheck; % size mismatch
for j = 1:m
    C_u_lo = [C_u_lo -sum(Dtau(:,j,:),3) + sum(reshape(DUcheck(:,j),r,s),2)<=0];
end
C = [Dtau >= 0,...
    C_bid_selection
    C_1agentpfield,...
    C_1fieldpagent,...
    C_agent_cap,...
    C_no_free_deli,...
    C_u_lo,...
    H_inventconstraints,...
    Atau_plus(:,:,q+1) == A_init,...
    H_invent(:,:,q+1) == C_init
];
ops = sdpsettings('solver','scip','verbose',0,'cachesolvers',1);
%ops.scip.maxnodes = 5e4;
%ops.scip.maxiter = 3e3;
%% global parameters
SetGlobal;
TTSUM = TSUMAN+TSUMMT; % simulation stop condition
%% get weather and season data
year_of_data = 1976;
station = 6;
yod = num2str(year_of_data);
weather_file =['ICWEATHR\NLD', num2str(station),'.', yod(2:4)];
%% season params
global DOYEM
% the day of the year on which crop emerges
DOYEM = 90;
% start time simulation
STTIME = 88;
% max season length
season_max_length = 300;
season_start_date(1:2) = monthday(STTIME); % season start: [M,D], e.g.: [3,4] equals March 4th
DOY = STTIME:season_max_length-1;          % day of year

%% weather data
[DTR,RAIN,TN,TX,~,VP,WN] = Read_Weatherfile(weather_file);
days = DOY(1):DOY(end);
DTR = DTR(days); RAIN = RAIN(days); TN = TN(days); TX = TX(days); VP = VP(days); WN = WN(days);
% DTR is daily total irradiation in MJ/m^2
% RAIN=max(0,RAIN+dRAIN(1:242)');
DAVTMP  = (TN+TX)/2;           % daily average temperature
DTEFF   = max(0,DAVTMP-TBASE); % daily effective temperature
%% initial conditions
[SLA_r,FLVT_r,NBALAN_r,TRANRF_r,TAGBM_r,RNSOIL_r,NTAC_r,NTAG_r,PEVAP_r,WATBAL_r,RWRT_r,RWST_r,LUECAL_r,RWSO_r,NDEMTO_r,RWLVG_r, FRTWET_r,FSTT_r,FSOT_r,FRACT_r,TSUM_r,NRF_r ...
    ,TEVAP_r,TTRAN_r,TRAIN_r,TRUNOF_r,TIRRIG_r,TDRAIN_r,TEXPLO_r,CUMPAR_r,GTSUM_r,WDRT_r,CBALAN_r,NLOSSR_r,NLOSSL_r,FERTNS_r,FERTN_r,TNSOIL_r,NUPTT_r,LAI_r,WLV_r,WLVG_r,WLVD_r,WSO_r,WST_r,WRT_r,ROOTD_r,WA_r ...
    ,WC_r,WCCR_r,DVS_r,NNI_r,ANLV_r,ANST_r,ANRT_r,ANSO_r,NUPTR_r]  =  deal(zeros(m,season_max_length));
% include random variables in intials.
SetInit;
%% get expected weather data
load('Expected_weather.mat');
load('Weather_mod.mat')
expected_weather = 0; % put this to zero to use old without expected weather.
DTEFF_E = zeros(length(DTEFF),1);
RAIN_E = zeros(length(RAIN),1);
DTR_E = zeros(length(DTR),1);
VP_E = zeros(length(VP),1);
WN_E = zeros(length(WN),1);
DAVTMP_E = zeros(length(DAVTMP),1);
n_w = 7;
%% Simulate
tic;
timer = toc;
dRAIN = zeros(1,season_max_length);
totU = zeros(2,season_max_length);
A = cell(season_max_length,1);
for t = 1:season_max_length-1
    % initialize bidding matrices
    G = zeros(s,m);
    L = zeros(r*s,m);
    %% prediction per client with or without expected weather
    if(expected_weather == 1)
        VP_E(t:t+n_w) = VP(t:t+n_w); VP_E(t+n_w+1:length(VP_E)) = mean_VP(t+n_w+1:length(VP_E));
        DTR_E(t:t+n_w) = DTR(t:t+n_w); DTR_E(t+n_w+1:length(DTR_E)) = mean_DTR(t+n_w+1:length(DTR_E));
        WN_E(t:t+n_w) = WN(t:t+n_w); WN_E(t+n_w+1:length(WN_E)) = mean_WN(t+n_w+1:length(WN_E));
        DAVTMP_E(t:t+n_w) = DAVTMP(t:t+n_w); DAVTMP_E(t+n_w+1:length(DAVTMP_E)) = mean_DAVTMP(t+n_w+1:length(DAVTMP_E));
        DTEFF_E = max(0,DAVTMP-TBASE);
        RAIN_E(t:t+n_w) = RAIN(t:t+n_w); RAIN_E(t+n_w+1:length(RAIN_E)) = mean_RAIN(t+n_w+1:length(RAIN_E));
    else
        DTEFF_E = DTEFF;
        VP_E = VP;
        DTR_E = DTR;
        WN_E = WN;
        DAVTMP_E = DAVTMP;
        RAIN_E = RAIN;
    end
    for  j = 1:m
        %% prediction per agent per scenario
        if TSUM_r(j,t) < TTSUM && DVS_r(j,t) < 2.01
            for k = 1:s
                [SLA,FLVT,NBALAN,TRANRF,TAGBM,RNSOIL,NTAC,NTAG,PEVAP,WATBAL,RWRT,RWST,LUECAL,RWSO,NDEMTO,RWLVG, ...
                    FRTWET,FSTT,FSOT,FRACT,TSUM,NRF,TEVAP,TTRAN,TRAIN,TRUNOF,TIRRIG,TDRAIN,TEXPLO,CUMPAR,GTSUM, ...
                    WDRT,CBALAN,NLOSSR,NLOSSL,FERTNS,FERTN,TNSOIL,NUPTT,LAI,WLVG,WLVD,WSO,WST,WRT,ROOTD,WA, ...
                    WC,WCCR,DVS,NNI,ANLV,ANST,ANRT,ANSO,NUPTR]  =  deal(zeros(1,season_max_length));
                % set initial conditions
                LAI(1)     = LAI_r(j,t);
                TSUM(1)    = TSUM_r(j,t);
                WC(1)      = WC_r(j,t);
                ROOTD(1)   = ROOTD_r(j,t);
                TEXPLO(1)  = TEXPLO_r(j,t);
                WA(1)      = WA_r(j,t);
                TEVAP(1)   = TEVAP_r(j,t);
                TTRAN(1)   = TTRAN_r(j,t);
                TRUNOF(1)  = TRUNOF_r(j,t);
                TIRRIG(1)  = TIRRIG_r(j,t);
                TDRAIN(1)  = TDRAIN_r(j,t);
                TRAIN(1)   = TRAIN_r(j,t);
                WLVG(1)    = WLVG_r(j,t);
                WST(1)     = WST_r(j,t);
                ANLV(1)    = ANLV_r(j,t);
                ANST(1)    = ANST_r(j,t);
                ANRT(1)    = ANRT_r(j,t);
                ANSO(1)    = ANSO_r(j,t);
                WRT(1)     = WRT_r(j,t);
                WSO(1)     = WSO_r(j,t);
                TNSOIL(1)  = TNSOIL_r(j,t);
                NUPTT(1)   = NUPTT_r(j,t);
                NLOSSL(1)  = NLOSSL_r(j,t);
                NLOSSR(1)  = NLOSSR_r(j,t);
                WLVD(1)    = WLVD_r(j,t);
                WDRT(1)    = WDRT_r(j,t);
                GTSUM(1)   = GTSUM_r(j,t);
                CUMPAR(1)  = CUMPAR_r(j,t);
                DVS1 = DVS1_r; DVS2 = DVS2_r; DSLR = DSLR_r;
                % next day after t
                i = 1;
                l = t;
                [TSUM(i+1),ROOTD(i+1),WA(i+1),WC(i+1),WCCR(i+1),TEXPLO(i+1),TEVAP(i+1),TTRAN(i+1),TRUNOF(i+1),TIRRIG(i+1),TDRAIN(i+1),TRAIN(i+1),DVS(i+1),NNI(i+1),SLA(i+1),LAI(i+1),...
                    NDEMTO(i+1),TNSOIL(i+1),NUPTT(i+1),ANLV(i+1),ANST(i+1),ANRT(i+1),ANSO(i+1),NLOSSL(i+1),NLOSSR(i+1),WLVG(i+1),WLVD(i+1),WST(i+1),WSO(i+1),WRT(i+1),TAGBM(i+1),...
                    NTAC(i+1),LUECAL(i+1),CUMPAR(i+1), GTSUM(i+1),WDRT(i+1),DSLR(j),DVS1(j),DVS2(j)]...
                    = LINTUL3_onestage(DOY(l),LAI(i),TSUM(i),DTEFF_E(l),WC(i),ROOTD(i),RAIN_E(l),DAVTMP_E(l),VP_E(l),DTR_E(l),WN_E(l),TEXPLO(i),WA(i),TEVAP(i),TTRAN(i),TRUNOF(i),TIRRIG(i),TDRAIN(i),TRAIN(i),WLVG(i),WST(i)...
                    ,ANLV(i),ANST(i),ANRT(i),ANSO(i),WRT(i),WSO(i),TNSOIL(i),NUPTT(i),NLOSSL(i),NLOSSR(i),WLVD(i),WDRT(i),GTSUM(i),CUMPAR(i),DSLR(j),DVS1(j),DVS2(j),U(1,k),0,U(2,k),0);
                i = i+1; l = l+1;
                % prediction of future days
                while l < length(DOY)-1 && TSUM(i) < TTSUM && DVS(i) < 2.01
                    % State update
                    [TSUM(i+1),ROOTD(i+1),WA(i+1),WC(i+1),WCCR(i+1),TEXPLO(i+1),TEVAP(i+1),TTRAN(i+1),TRUNOF(i+1),TIRRIG(i+1),TDRAIN(i+1),TRAIN(i+1),DVS(i+1),NNI(i+1),SLA(i+1),LAI(i+1),...
                        NDEMTO(i+1),TNSOIL(i+1),NUPTT(i+1),ANLV(i+1),ANST(i+1),ANRT(i+1),ANSO(i+1),NLOSSL(i+1),NLOSSR(i+1),WLVG(i+1),WLVD(i+1),WST(i+1),WSO(i+1),WRT(i+1),TAGBM(i+1),...
                        NTAC(i+1),LUECAL(i+1),CUMPAR(i+1), GTSUM(i+1),WDRT(i+1),DSLR(j),DVS1(j),DVS2(j)]...
                        = LINTUL3_onestage(DOY(l),LAI(i),TSUM(i),DTEFF_E(l),WC(i),ROOTD(i),RAIN_E(l),DAVTMP_E(l),VP_E(l),DTR_E(l),WN_E(l),TEXPLO(i),WA(i),TEVAP(i),TTRAN(i),TRUNOF(i),TIRRIG(i),TDRAIN(i),TRAIN(i),WLVG(i),WST(i)...
                        ,ANLV(i),ANST(i),ANRT(i),ANSO(i),WRT(i),WSO(i),TNSOIL(i),NUPTT(i),NLOSSL(i),NLOSSR(i),WLVD(i),WDRT(i),GTSUM(i),CUMPAR(i),DSLR(j),DVS1(j),DVS2(j),0,0,0,0);
                    % increment day
                    i = i+1;
                    l = l+1;
                end % end prediction
                G(k,j) = WSO(i-1)-WSO_r(j,t);
            end % end bids
        end
    end  % end clients
    %% Selection
    j_g = kron(delt, pi)'*G;  %kron(ones(s,1),pi)'*(G.*kron(delt,ones(p,1)))*ones(m,1);
    j_l = kron(delt,rho)'*L;  %kron(ones(s,1),rho)'*(L.*kron(delt,ones(r,1)))*ones(m,1);    
    for z = 1:q
        if z==1
            j_m1 = trace(Atau_plus(:,:,q+1)*E*Atau_plus(:,:,z)');
        else
            j_m2 = trace(Atau_plus(:,:,z-1)*E*Atau_plus(:,:,z)');
        end
    end
    j_q = kron(ones(m,1),rho')*sum(Dtau,3)+ kron(ones(m,1),sigma')*sum(Atau_plus(:,1:m,:),3)+j_m1+j_m2;% rho'*sum(sum(Dtau,2),3) + sigma'*sum(sum(Atau,2),3);
    J = trace(j_g - j_l - j_q);
    optimize(C, -J ,ops);    
    D = round(value(delt));
    A{t} = round(value(Atau_plus)); 
    H{t} = round(value(H_invent),5); % use round here because scip value can be ~~ 0 example: -0.00000001, if this initial condition is used the 
    % the solution fails 
    C_input{t} = value(C_in);
    Dtau_out{t} = round(value(Dtau));
    C(end) =[]; % delete previous initial inventory
    C(end) =[]; % delete previous initial position agents
    C = [C, Atau_plus(:,:,q+1) == A{t}(:,:,q)]; % add initial condition agents
    C = [C, H_invent(:,:,q+1) == H{t}(:,:,q)];
  
    %% update
    if(expected_weather == 1)
        dRAIN = std_RAIN(t)*(randn-1);
        % some variables are correlated so
        vec = [std_DTR(k) zeros(1,3); 0 std_VP(k) 0 0; 0 0 std_WN(k) 0; 0 0 0 std_DAVTMP(k)]*B_weather_mod{t+STTIME}*rand(4,1);
        dDTR = vec(1); dVP = vec(2); dWN = vec(3); dDAVTMP = vec(4);
    else
        dRAIN = 0;
        dDTEFF = 0;
        dWN  = 0;
        dVP = 0;
        dDAVTMP = 0;
        dDTR = 0;
    end
    totU(:,t) = sum(U*D,2);
    for j = 1:m
        if TSUM_r(j,t) < TTSUM && DVS_r(j,t) < 2.01
            U_temp = U*D(:,j); U_irrig(j,t) = U_temp(1); U_fert(j,t) = U_temp(2);            
            [TSUM_r(j,t+1),ROOTD_r(j,t+1),WA_r(j,t+1),WC_r(j,t+1),WCCR_r(j,t+1),TEXPLO_r(j,t+1),TEVAP_r(j,t+1),TTRAN_r(j,t+1),TRUNOF_r(j,t+1),TIRRIG_r(j,t+1),TDRAIN_r(j,t+1),TRAIN_r(j,t+1),DVS_r(j,t+1),NNI_r(j,t+1),SLA_r(j,t+1)...
                ,LAI_r(j,t+1),NDEMTO_r(j,t+1),TNSOIL_r(j,t+1),NUPTT_r(j,t+1),ANLV_r(j,t+1),ANST_r(j,t+1),ANRT_r(j,t+1),ANSO_r(j,t+1),NLOSSL_r(j,t+1),NLOSSR_r(j,t+1),WLVG_r(j,t+1),WLVD_r(j,t+1),WST_r(j,t+1),WSO_r(j,t+1)...
                ,WRT_r(j,t+1),TAGBM_r(j,t+1),NTAC_r(j,t+1),LUECAL_r(j,t+1),CUMPAR_r(j,t+1), GTSUM_r(j,t+1),WDRT_r(j,t+1),DSLR_r(j),DVS1_r(j),DVS2_r(j)]...
                = LINTUL3_onestage(DOY(t),LAI_r(j,t),TSUM_r(j,t),max(0,DAVTMP(t)+dDAVTMP-TBASE),WC_r(j,t),ROOTD_r(j,t),max(0,dRAIN+RAIN(t)),DAVTMP(t)+dDAVTMP,max(0,VP(t)+dVP),max(0,DTR(t)+dDTR),max(0,WN(t)+dWN),TEXPLO_r(j,t),WA_r(j,t),TEVAP_r(j,t),TTRAN_r(j,t),TRUNOF_r(j,t)...
                ,TIRRIG_r(j,t),TDRAIN_r(j,t),TRAIN_r(j,t),WLVG_r(j,t),WST_r(j,t),ANLV_r(j,t),ANST_r(j,t),ANRT_r(j,t),ANSO_r(j,t),WRT_r(j,t),WSO_r(j,t),TNSOIL_r(j,t),NUPTT_r(j,t),NLOSSL_r(j,t),NLOSSR_r(j,t)...
                ,WLVD_r(j,t),WDRT_r(j,t),GTSUM_r(j,t),CUMPAR_r(j,t),DSLR_r(j),DVS1_r(j),DVS2_r(j),U_irrig(j,t),0,U_fert(j,t),0);
        else
            TSUM_r(j,t+1) = TSUM_r(j,t);
            ROOTD_r(j,t+1)= ROOTD_r(j,t);
            WA_r(j,t+1) = WA_r(j,t);
            WC_r(j,t+1) = WC_r(j,t);
            WCCR_r(j,t+1) = WCCR_r(j,t);
            TEXPLO_r(j,t+1) = TEXPLO_r(j,t);
            TEVAP_r(j,t+1) = TEVAP_r(j,t);
            TTRAN_r(j,t+1) = TTRAN_r(j,t);
            TRUNOF_r(j,t+1) = TRUNOF_r(j,t);
            TIRRIG_r(j,t+1) = TIRRIG_r(j,t);
            TDRAIN_r(j,t+1) = TDRAIN_r(j,t);
            TRAIN_r(j,t+1) = TRAIN_r(j,t);
            DVS_r(j,t+1) = DVS_r(j,t);
            NNI_r(j,t+1) = NNI_r(j,t);
            SLA_r(j,t+1) = SLA_r(j,t);
            LAI_r(j,t+1) = LAI_r(j,t);
            NDEMTO_r(j,t+1) = NDEMTO_r(j,t);
            TNSOIL_r(j,t+1) = TNSOIL_r(j,t);
            NUPTT_r(j,t+1) = NUPTT_r(j,t);
            ANLV_r(j,t+1) = ANLV_r(j,t);
            ANST_r(j,t+1) = ANST_r(j,t);
            ANRT_r(j,t+1) = ANSO_r(j,t);
            ANSO_r(j,t+1) = ANSO_r(j,t);
            NLOSSL_r(j,t+1) = NLOSSL_r(j,t);
            NLOSSR_r(j,t+1) = NLOSSR_r(j,t);
            WLVG_r(j,t+1) = WLVG_r(j,t);
            WLVD_r(j,t+1) = WLVD_r(j,t);
            WST_r(j,t+1) = WST_r(j,t);
            WSO_r(j,t+1) = WSO_r(j,t);
            WRT_r(j,t+1) = WRT_r(j,t);
            TAGBM_r(j,t+1) = TAGBM_r(j,t);
            NTAC_r(j,t+1) = NTAC_r(j,t);
            LUECAL_r(j,t+1) = LUECAL_r(j,t);
            CUMPAR_r(j,t+1) = CUMPAR_r(j,t);
            GTSUM_r(j,t+1) = GTSUM_r(j,t);
            WDRT_r(j,t+1) = WDRT_r(j,t);
        end
    end
    dt = toc - timer;
    timer = timer+dt;
    e_t = 0.5*dt*(i+1);
    fprintf('%u/%u = %0.2f%% - Runtime %0.1f [s] - Expected duration %u [m] %0.0f [s]\n',t,season_max_length,100*t/season_max_length,timer,floor(e_t/60),rem(e_t,60));
    if sum(TSUM_r(:,t)>=TTSUM) == m || sum(DVS_r(:,t)>=2.01) == m
        break
    end
end
% rescale variables add lintul 3 ones
DOY     = DOY(:,1:t);
TSUM    = TSUM_r(:,1:t);
LAI     = LAI_r(:,1:t);
WLV     = WLV_r(:,1:t);
WLVG    = WLVG_r(:,1:t);
WLVD    = WLVD_r(:,1:t);
WSO     = WSO_r(:,1:t);
WST     = WST_r(:,1:t);
WRT     = WRT_r(:,1:t);
ROOTD   = ROOTD_r(:,1:t);
WA      = WA_r(:,1:t);
WC      = WC_r(:,1:t);
WSOTHA  = WSO/100;
prof = pi*sum(WSO(:,end)) - sum(rho.*sum(totU,2))

%% save results
Str = strcat('results/DHMARA_refill_t_m', num2str(m),'_p', num2str(p),'_q',num2str(q),'_r',num2str(r),'_n',num2str(n));
CheckStr = strcat(Str,'.mat');
if(exist(CheckStr,'file') == 2)
    prompt = 'file already exists, save anyway y/n?';
    x = input(prompt,'s');
    if(x == 'y')
        disp('saved')
        save(Str)
    end
else
    save(Str)
end
%% plot
figure(1)
subplot(413)
plotyear(DOY(1:length(U_fert)),sum(U_fert));
ylabel(['$\Sigma_{\j \in M} \quad fertilization$'],'interpreter','latex')
subplot(414)
plotyear(DOY(1:length(U_irrig)),sum(U_irrig));
xlabel('Day of year')
ylabel(['$\Sigma_{\j \in M} \quad irrigation$'],'interpreter','latex')
subplot(4,1,[1 2])
title_str = strcat("Number of agents = ", num2str(n),", number of client = ", num2str(m));
plotyear(DOY,WSOTHA,LAI);
title(title_str)
grid on

%% save for gui
count = 1;
while(true)
    GuiStr = ['results/gui_refill' num2str(count) '.mat'];
    if(exist(GuiStr,'file') == 2)
        % dont save
        count = count+1;
    else
        save(GuiStr,'DOY', 'DTR', 'WN','DAVTMP','DTEFF','RAIN','U_irrig','U_fert','TSUM_r','ROOTD_r','WA_r','WC_r','WCCR_r','TEXPLO_r', 'TEVAP_r','TTRAN_r','TRUNOF_r','TIRRIG_r'...
            ,'TDRAIN_r','TRAIN_r','DVS_r', 'NNI_r','SLA_r','LAI_r','NDEMTO_r','TNSOIL_r','NUPTT_r','ANLV_r','ANST_r','ANRT_r','ANSO_r','NLOSSL_r','NLOSSR_r','WLVG_r','WST_r','WSO_r','WRT_r','TAGBM_r', ...
            'NTAC_r','LUECAL_r','CUMPAR_r','GTSUM_r','WDRT_r','Cmax','q','p','A');
        break;
    end
end







