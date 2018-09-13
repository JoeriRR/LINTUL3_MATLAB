clear variables; close all;
%% plot
% read res.dat LINTUL 3 from fst simulator (wageningen university)
[time,~,DVS_fst,TSUM_fst,TAGBM_fst,WST_fst,WLVG_fst,WLVD_fst,WSO_fst,LAI_fst,NTAC_fst,WRT_fst,GTSUM_fst, CBALAN_fst, TRANRF_fst, NNI_fst, SLA_fst, FRACT_fst,FRTWET_fst,FLVT_fst,FSTT_fst,FSOT_fst,RWLVG_fst,RWST_fst,RWRT_fst,RWSO_fst,CUMPAR_fst,LUECAL_fst,NUPTT_fst,TTRAN_fst,TEVAP_fst,PEVAP_fst,NBALAN_fst,WATBAL_fst,NUPTR_fst,TNSOIL_fst,NDEMTO_fst,RNSOIL_fst,FERTN_fst,FERTNS_fst,WA_fst,TIRRIG_fst,TRAIN_fst,TEXPLO_fst,TRUNOF_fst,TDRAIN] = importfile('res.dat'); 
% put all variables in a matrix
A = [DVS_fst,TSUM_fst,TAGBM_fst,WST_fst,WLVG_fst,WLVD_fst,WSO_fst,LAI_fst,NTAC_fst,WRT_fst,GTSUM_fst, CBALAN_fst, TRANRF_fst, NNI_fst, SLA_fst, FRACT_fst,FRTWET_fst,FLVT_fst,FSTT_fst,FSOT_fst,RWLVG_fst,RWST_fst,RWRT_fst,RWSO_fst,CUMPAR_fst,LUECAL_fst,NUPTT_fst,TTRAN_fst,TEVAP_fst,PEVAP_fst,NBALAN_fst,WATBAL_fst,NUPTR_fst,TNSOIL_fst,NDEMTO_fst,RNSOIL_fst,FERTN_fst,FERTNS_fst,WA_fst,TIRRIG_fst,TRAIN_fst,TEXPLO_fst,TRUNOF_fst,TDRAIN];
% name vector
names = ["DVS","TSUM","TAGBM","WST","WLVG","WLVD","WSO","LAI","NTAC","WRT","GTSUM"," CBALAN"," TRANRF"," NNI"," SLA"," FRACT","FRTWET","FLVT","FSTT","FSOT","RWLVG","RWST","RWRT","RWSO","CUMPAR","LUECAL","NUPTT","TTRAN","TEVAP","PEVAP","NBALAN","WATBAL","NUPTR","TNSOIL","NDEMTO","RNSOIL","FERTN","FERTNS","WA","TIRRIG","TRAIN","TEXPLO","TRUNOF","TDRAIN"];
figcount = 0;
% plot
vec = 1:240;
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
% load matlab LINTUL 3 variables
load('Z.mat')
% put all variables in matrix
B = [DVS,TSUM,TAGBM,WST,WLVG,WLVD,WSO,LAI,NTAC,WRT,GTSUM, CBALAN, TRANRF, NNI, SLA, FRACT,FRTWET,FLVT,FSTT,FSOT,RWLVG,RWST,RWRT,RWSO,CUMPAR,LUECAL,NUPTT,TTRAN,TEVAP,PEVAP,NBALAN,WATBAL,NUPTR,TNSOIL,NDEMTO,RNSOIL,FERTN,FERTNS,WA,TIRRIG,TRAIN,TEXPLO,TRUNOF,TDRAIN];
% plot over previous figure to compare
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
clc;
% prints the 2-norm of all saved variables.
for k = 1:length(names)
   temp = norm(A(vec,k)-B(vec,k),2);
   NormStr = strcat("Var = ",names(k)," 2-norm(Var_FST-Var_MATLAB)) =", num2str(temp));
   disp(NormStr)
end


