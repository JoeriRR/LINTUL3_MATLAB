%% load
clear variables; close all; clc;
% read res.dat LINTUL 3 from fst simulator (wageningen university)
[time,~,DVS_fst,TSUM_fst,TAGBM_fst,WST_fst,WLVG_fst,WLVD_fst,WSO_fst,LAI_fst,NTAC_fst,WRT_fst,GTSUM_fst, CBALAN_fst, TRANRF_fst, NNI_fst, SLA_fst, FRACT_fst,FRTWET_fst,FLVT_fst,FSTT_fst,FSOT_fst,RWLVG_fst,RWST_fst,RWRT_fst,RWSO_fst,CUMPAR_fst,LUECAL_fst,NUPTT_fst,TTRAN_fst,TEVAP_fst,PEVAP_fst,NBALAN_fst,WATBAL_fst,NUPTR_fst,TNSOIL_fst,NDEMTO_fst,RNSOIL_fst,FERTN_fst,FERTNS_fst,WA_fst,TIRRIG_fst,TRAIN_fst,TEXPLO_fst,TRUNOF_fst,TDRAIN_fst] = importfile('res.dat');
% put all variables in a matrix
A = [DVS_fst,TSUM_fst,TAGBM_fst,WST_fst,WLVG_fst,WLVD_fst,WSO_fst,LAI_fst,NTAC_fst,WRT_fst,GTSUM_fst, CBALAN_fst, TRANRF_fst, NNI_fst, SLA_fst, FRACT_fst,FRTWET_fst,FLVT_fst,FSTT_fst,FSOT_fst,RWLVG_fst,RWST_fst,RWRT_fst,RWSO_fst,CUMPAR_fst,LUECAL_fst,NUPTT_fst,TTRAN_fst,TEVAP_fst,PEVAP_fst,NBALAN_fst,WATBAL_fst,NUPTR_fst,TNSOIL_fst,NDEMTO_fst,RNSOIL_fst,FERTN_fst,FERTNS_fst,WA_fst,TIRRIG_fst,TRAIN_fst,TEXPLO_fst,TRUNOF_fst,TDRAIN_fst];
% load matlab LINTUL 3 variables
load('results\Z.mat')
% put all variables in matrix
B = [DVS,TSUM,TAGBM,WST,WLVG,WLVD,WSO,LAI,NTAC,WRT,GTSUM, CBALAN, TRANRF, NNI, SLA, FRACT,FRTWET,FLVT,FSTT,FSOT,RWLVG,RWST,RWRT,RWSO,CUMPAR,LUECAL,NUPTT,TTRAN,TEVAP,PEVAP,NBALAN,WATBAL,NUPTR,TNSOIL,NDEMTO,RNSOIL,FERTN,FERTNS,WA,TIRRIG,TRAIN,TEXPLO,TRUNOF,TDRAIN];
% name vector
names = ["DVS","TSUM","TAGBM","WST","WLVG","WLVD","WSO","LAI","NTAC","WRT","GTSUM","CBALAN","TRANRF","NNI","SLA","FRACT","FRTWET","FLVT","FSTT","FSOT","RWLVG","RWST","RWRT","RWSO","CUMPAR","LUECAL","NUPTT","TTRAN","TEVAP","PEVAP","NBALAN","WATBAL","NUPTR","TNSOIL","NDEMTO","RNSOIL","FERTN","FERTNS","WA","TIRRIG","TRAIN","TEXPLO","TRUNOF","TDRAIN"];
figcount = 0;
vec = 1:240;
% plot all variables
for k = 1:length(names)
    if(mod(k,9)==1)
        figcount = figcount+1;
        figure();
    end
    figure(figcount)
    subplot(3,3,k-((figcount-1)*9))
    plot(vec,A(:,k)); hold on;
    ylabel(names(k));
    xlabel('day')
end
figcount = 0;
for k = 1:length(names)
    if(mod(k,9)==1)
        legend('FST', 'MATLAB')
        figcount = figcount+1;
    end
    figure(figcount)
    subplot(3,3,k-((figcount-1)*9))
    plot(vec,B(vec,k));
end

%% compare variables
clc; day = zeros(size(A,2),1);
% prints the Mean squared error
ndigits = 10;
for k = 1:length(names)
    temp = norm(A(vec,k)-B(vec,k),2)/length(vec); % MSE
    [temp2,I]= max(abs(A(vec,k)-B(vec,k))); % max daily error
    if(max(A(vec,k))~= 0)
        temp3 = temp2/max(abs(A(vec,k)))*100; % max percentage error
    end
    day(k) = I;
    if(I == 1)
        I = [];
    end
    NormStr = pad(strcat(pad(names(k),10),"MSE = ", num2str(temp)),28);
    MaxStr = pad(strcat(pad(strcat("MAX error daily = ",num2str(temp2)),30),"At day ", num2str(I)),40);
    PercStr = strcat("  Max % & ", num2str(temp3));
    PrintStr = strcat(NormStr, MaxStr, PercStr);  
    digits(ndigits)
    LatexTab = strcat(names(k),'&',char(vpa(round(temp,ndigits))),'&', char(vpa(round(temp2,ndigits))),'&',char(vpa(round(temp3,ndigits))),'&',num2str(I),'\\');
    %disp(LatexTab)
    disp(PrintStr)
end

%% important variable plots
close all
figure(1);clf
%plotyear(time,WSO(vec)/100,LAI(vec));
plot(time,WSO(vec)/100, 'Linewidth',2)
xlabel('Day of year')
hold on
plot(time,LAI(vec),'linewidth',2)
plot(time,WSO_fst(vec)/100,'--','LineWidth',2)
plot(time,LAI_fst(vec),'--','LineWidth',2);
xlim([1 365]);grid on;
    set(gca,'XTick',[1,32,60,91,121,152,182,213,244,274,305,335])
    set(gca,'XTickLabel',{'            Jan','            Feb','            Mar',...
        '            Apr','            May','            Jun','            Jul',...
        '            Aug','            Sep','            Oct','            Nov',...
        '            Dec'})
%plotyear(time,WSO_fst/100,LAI_fst)
legend('MATLAB WSOTHA','MATLAB LAI', 'FST WSOTHA', 'FST LAI')
grid on






