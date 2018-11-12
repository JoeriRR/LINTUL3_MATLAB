clear variables; clc;
% script to calculate table 4
addpath('Init_scripts');
E = GetE_gui_grid(20,1);
load('results/DHMARA_m20_p1_q2_r2_n2.mat');
dist_old = 0;
dist_new = 0;
A_init_pos = [eye(n) zeros(n,m-n)]; % assume same initial condition
% table report: files DHMARA_spatial_m20_p1_q2_r2_n2
% table report: files DHMARA_spatial_largeE_m20_p1_q2_r2_n2
% table report: files DHMARA_spatial_verylargeE_m20_p1_q2_r2_n2
for i = 1:length(A_temp)
    if(isempty(A_temp{i}))
        break;
    end
    for j = 1:q
        % assume idle agents do not move
        for x = 1:n
            if(isequal(A_temp{i}(x,:,j),zeros(1,m)) == 1)
                A_temp{i}(x,:,j) = A_init_pos(x,:);
            end            
        end
        if(j == 1)
            dist_old = dist_old+trace(A_init_pos*E*A_temp{i}(:,:,j)');
        else                       
            dist_old = dist_old+trace(A_temp{i}(:,:,j-1)*E*A_temp{i}(:,:,j)');
        end
        % update initial condition
        for x = 1:n
            if(isequal(A_temp{i}(x,:,j),zeros(1,m)) == 0)
                A_init_pos(x,:) = A_temp{i}(x,:,j);
            end
        end
    end
end
dist_old
WSO_old = sum(WSO(:,end))
n_irrig_old = length(find(U_irrig>0))
n_fert_old = length(find(U_fert>0))
load('results/DHMARA_spatial_m20_p1_q2_r2_n2.mat');
for i = 1:length(A)
    if(isempty(A{i}))
        break;
    end
    for j = 1:q
        if(j==1)
            dist_new = dist_new + trace(A{i}(:,:,3)*E*A{i}(:,:,j)');
        else
            dist_new = dist_new + trace(A{i}(:,:,j-1)*E*A{i}(:,:,j)');
        end
    end
end
dist_new
WSO_new = sum(WSO(:,end))
n_irrig_new = length(find(U_irrig>0))
n_irrig_fert = length(find(U_fert>0))