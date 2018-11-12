clear variables; close all; clc;
load('results/DHMARA_refill_m5_p1_q3_r2_n2.mat');
% figure 7 from report, change loadname and nn to 20 to get a better plot.
figure(1)
nn = 130;
var_in1 = zeros(nn*q,n);
var_in11 = var_in1;
var_in111 = var_in11;
index = 1;
for t = 1:nn
    for tau = 1:q
        for i = 1:n
            tempvar = find(round(H{t}(:,i,tau),2)); %  index inventory agent i
            if(isempty(tempvar))
                var_in1(index,i) = 0;
            else
                var_in1(index,i) = H{t}(tempvar,i,tau); % inventory agent i
            end
            var_in11(index,i) = find(A{t}(i,:,tau)); % allocation agent i
            Z = Dtau_out{t}(:,:,tau)*A{t}(i,1:m,tau)';
            P = find(Z);
            if(isempty(P))
                var_in111(index,i) =0 ;
            else               
                var_in111(index,i) = Z(P);
            end           
        end
        index = index +1;
    end
end
xtickss = 1:q:nn*q;
for i = 1:length(xtickss)
    labels{i} = num2str(i);
end
for i = 1:n
    figure(i);
    subplot(311)
    stairs(1:(index-1),var_in1(:,i),'linewidth',2);
    hold on
    aaaa = find(Cmax(:,i));
    plot(1:(index-1),Cmax(aaaa,i)*ones(index-1,1),'Linestyle','- -','linewidth',2);
    title(['agent ' num2str(i) ', m = ' num2str(m) ', q = ' num2str(q)])
    legend('H', 'max capacity')
    ylim([-2 Cmax(aaaa,i)*1.2])
    xlim([1 nn*q])
    set(gca,'xtick', xtickss)
    set(gca,'xticklabel',labels)
    ylabel('inventory')
    subplot(313)
    stairs(1:(index-1),var_in11(:,i),'linewidth',2);
    ylabel('allocation')
    ylim([0 m+3]);
    xlim([1 nn*q])
    xlabel('t [days]')
    set(gca,'xtick', xtickss)
    set(gca,'xticklabel',labels)
    subplot(312)
    stairs(1:(index-1),var_in111(:,i),'linewidth',2);
    ylabel('delivery');
    ylim([-1 12])
    xlim([1 nn*q])
    set(gca,'xtick', xtickss)
    set(gca,'xticklabel',labels)
    hold off;
end

