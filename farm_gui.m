function farm_gui
%% create figure, buttons and load simulated results
close all; clc;
data = load('results/gui.mat');
names = {'DTR', 'WN','DAVTMP','DTEFF','RAIN','U_irrig','U_fert','TSUM_r','ROOTD_r','WA_r','WC_r','WCCR_r','TEXPLO_r', 'TEVAP_r','TTRAN_r','TRUNOF_r','TIRRIG_r'...
    ,'TDRAIN_r','TRAIN_r','DVS_r', 'NNI_r','SLA_r','LAI_r','NDEMTO_r','TNSOIL_r','NUPTT_r','ANLV_r','ANST_r','ANRT_r','ANSO_r','NLOSSL_r','NLOSSR_r','WLVG_r','WST_r','WSO_r','WRT_r','TAGBM_r', ...
    'NTAC_r','LUECAL_r','CUMPAR_r','GTSUM_r','WDRT_r'};
%  Create and then hide the GUI as it is being constructed, figure width:
farmlength = 10;
agentlength = 2;
SwitchCase = 0;
% asssign varibles ie qmax,m,n,radius etc. create function for this
m = size(data.LAI_r,1); radius = 0.2;
Cmax = [5 5 0; 0 0 10]; % agent max capacity matrix four agents that can deliver water and 1 agents that can deliver fertilizer
[r,n] = size(Cmax);
for b=1:m
    fieldnames{b} = num2str(b);
end
fieldnames{b+1} = 'all';
[end_date] = length(data.U_irrig); q_sim = 2; % need to add U and delta matrix to assign the agents!
%% buttons and text
%position [left bottom width heigth]!
f = figure('Visible','off','Position',[100,200,900,600]);
set(f,'MenuBar','none')
set(f,'ToolBar','none')
figend = 670; buttonwidth = 100; figtop = 525;
%hsave =
hsave = uicontrol('Style','pushbutton','String','save','Position',[figend+110,figtop-250,buttonwidth,25],'Callback', {@savebutton_Callback});
hplay = uicontrol('Style','togglebutton','String','play','Position',[figend,figtop,buttonwidth,25],'Callback',{@playbutton_Callback});
hplot = uicontrol('Style','togglebutton','String','plot','Position',[figend+110,figtop,buttonwidth,25],'Callback',{@plotbutton_Callback});
hdata = uicontrol('Style','popupmenu', 'String',names,'Position',[figend+10,figtop-60,buttonwidth-20,25],'Callback',{@varpopup_Callback});
fdata = uicontrol('Style','popupmenu', 'String',fieldnames,'Position',[figend+120,figtop-60,buttonwidth-20,25],'Callback',{@fieldpopup_Callback});
figoptions = uicontrol('Style','text','String','select variable:','Position',[figend, figtop-35, buttonwidth, 25],'Fontsize',9);
figfield = uicontrol('Style','text','String','select field:','Position',[figend+110, figtop-35, buttonwidth, 25],'Fontsize',9);

slider_day = uicontrol('Style','slider','Position',[80,20,550,25],'min',1,'max',end_date,'Callback',{@sliderbutton_Callback});
set(slider_day,'Value', 1);
StaticField = uicontrol('Style','text', 'String','field status:','Position', [figend, figtop-140, 100, 15]);
StaticField2 = uicontrol('Style','text', 'String','agent status:','Position', [figend+110, figtop-140, 100, 15]);
dayweathertext = uicontrol('Style','text','String','DOY','Position',[80, 550, 300,25], 'fontsize', 11);
fieldtext = uicontrol('Style','text','String','-','Position',[figend, figtop-200, 100, 50]);
agenttext = uicontrol('Style','text','String','-','Position',[figend+110, figtop-200, 100, 50]);
ha = axes('Units','pixels','Position',[80,100,550,450],'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
f.Visible = 'on';

%% initialise UI and make all buttons resizable
f.Units = 'normalized';
ha.Units = 'normalized';
hplay.Units = 'normalized';
fieldtext.Units = 'normalized';
agenttext.Units = 'normalized';
StaticField.Units = 'normalized';
StaticField2.Units = 'normalized';
dayweathertext.Units = 'Normalized';
hplot.Units = 'normalized';
fdata.Units = 'normalized';
figoptions.Units = 'normalized';
figfield.Units = 'normalized';
slider_day.Units = 'normalized';
hdata.Units = 'normalized';
hsave.Units = 'normalized';
range = [-agentlength farmlength 0 farmlength];
axis(range);
%% main
while(true)
    if(isvalid(f))
        switch SwitchCase
            case 0
                k = 1; % day
                q = 1; % subtime unit
                [center, fieldwidth,fieldheigth]= draw_Field(m,range,5/7);
                [circles,init_c_agent] = init_Agent(n,range);
                center_agent = init_c_agent; draw_Border(range);
                set_Day_Weather(data);
                id1 = iptaddcallback(f,'WindowButtonMotionFcn',@selecter);
                Current_x = data.DOY; Current_y = data.LAI_r; varname ='LAI';
                field = 1; last_entry = false; data_changed = true;
                disp('Initialise');
                SwitchCase = 1;
            case 1
                if(get(hplay,'Value'))
                    % increment day
                    disp('playing data')
                    if (k >= end_date)
                        k = end_date;
                    else
                        if(mod(q,q_sim) == 0)
                            q = 0;
                            k = k+1;
                            data_changed = true;
                        end
                        q = q+1;
                    end
                    set(slider_day,'value',k);
                else
                    k = round(get(slider_day,'value'));
                end
                if(get(hplot,'Value'))
                    if(last_entry == true)
                        iptremovecallback(f,'WindowButtonMotionFcn', id1);
                        grid on;
                    end
                    % check if data changed
                    if(data_changed == true)
                        cla;
                        disp('plotting data')
                        plotdata(0);
                    end
                    data_changed = false;
                    last_entry = false;
                else
                    % only do once
                    if(last_entry == false)
                        id1 = iptaddcallback(f,'WindowButtonMotionFcn',@selecter);
                        axis off;
                        cla;
                        legend('hide')
                        [center, fieldwidth,fieldheigth]= draw_Field(m,range,5/7);
                        [indecis] = assign_agent(Cmax);
                        [circles,center_agent] = move_agent(circles,indecis);
                        draw_Border(range);
                        axis(range);
                    else
                        
                        disp('moving agents')
                        [indecis] = assign_agent(Cmax);
                        [circles,center_agent] = move_agent(circles,indecis);
                        
                    end
                    last_entry = true;
                end
                set_Day_Weather(data);
                pause(0.5);
        end
    else
        return; % close while loop when exit button is pressed
    end
end
%% button functions
    function playbutton_Callback(source,eventdata)
        if(get(hplay,'Value'))
            set(hplay,'String', 'pause')
        else
            set(hplay,'String', 'play')
        end
    end
    function plotbutton_Callback(source,eventdata)
        if(get(hplot,'Value') == 1)
            set(hplot,'String', 'vis')
        else
            set(hplot,'String','plot')
        end
        data_changed = true;
    end
    function varpopup_Callback(hdata,source,eventdata)
        val = hdata.Value;
        Str = hdata.String;
        varname = Str{val};
        Current_y = data.(varname);
        data_changed =true;
    end
    function fieldpopup_Callback(fdata,source,eventdata)
        val = fdata.Value;
        Str = fdata.String;
        if(strcmp(Str{val},'all'))
            field = m+1;
        else
            field = str2double(Str{val});
        end
        data_changed = true;
    end
    function sliderbutton_Callback(hplay,source,eventdata)
        data_changed = true;
    end
    function savebutton_Callback(source,eventdata)
        if(get(hplot,'value'))
            plotdata(1);
        end
    end
%% textbox functions
    function set_Day_Weather(data)
        DOY = data.DOY;
        RAIN = data.RAIN;
        set(dayweathertext,'String',['day ' num2str(DOY(k)) ' subtime unit q ' num2str(q) ' rainfall ' num2str(RAIN(k)) '[mm]']);
    end
    function [current_pos] = selecter(source,eventdata)
        current_pos = get(gca,'CurrentPoint');
        for i = 1:size(center,1)
            if(current_pos(1,1)>=(center(i,1)-fieldwidth/2) && current_pos(1,1)<= (center(i,1)+fieldwidth/2) && current_pos(1,2)>=(center(i,2)-fieldheigth/2) && current_pos(1,2)<= (center(i,2)+fieldheigth/2))
                set(fieldtext,'String',['Field ',num2str(i), newline, ['LAI ' num2str(data.LAI_r(i,k))], newline, ['WSO ' num2str(data.WSO_r(i,k))]]);
                break;
            else
                %set(fieldtext,'String','');
            end
        end
        for i = 1:size(center_agent,1)
            r = sqrt((current_pos(1,1)-center_agent(i,1))^2+(current_pos(1,2)-center_agent(i,2))^2);
            if(r<=radius)
                set(agenttext,'String',['Agent ',num2str(i), newline, 'res 1', newline, 'res 2']);
                break;
            else
                %set(agenttext,'String','');
            end
        end
    end

%% draw functions
    function [center,fieldwidth,fieldheigth] = draw_Field(m,range,fillfact)
        % draw fields as rectangles
        figure(f);
        rows = round(sqrt(m));
        if rows^2 < m
            rows = rows+1;
        end
        if (mod(m,rows) == 0)
            columns = m/rows;
        else
            columns = rows;
        end
        fieldwidth = abs(0-range(2))*fillfact/columns;
        openwidth = (abs(0-range(2))-fieldwidth*columns)/(columns+1);
        fieldheigth = abs(range(3)-range(4))*fillfact/rows;
        openheigth = (abs(range(3)-range(4))-fieldheigth*rows)/(rows+1);
        index = 1; center = zeros(m,2);
        for i = 1:rows
            for j = 1:columns
                rectangle('Position',[openwidth*j+(j-1)*fieldwidth openheigth*i+(i-1)*fieldheigth fieldwidth fieldheigth]);
                center(index,1:2) = [(openwidth*j+(j-1)*fieldwidth+fieldwidth/2) openheigth*i+(i-1)*fieldheigth+fieldheigth/2];
                if index == m
                    return;
                end
                index = index+1;
            end
        end
    end
    function draw_Border(range)
        rectangle('Position', [range(1) range(3) abs(range(1)-range(2)) abs(range(3)-range(4))],'Linestyle','--');
    end
    function [circleobj,center] = init_Agent(n,range)
        index = 1;
        gap = (abs(range(3)-range(4))-radius*n*2)/(n+1);
        center = zeros(n,2);
        for i = 1:n
            center(index,1:2) = [range(1)/2, gap*i+(radius*(2*i-1))];
            index = index+1;
        end
        circleobj = viscircles(center,radius*ones(n,1));
    end
    function [circleobjout,cur_agent_pos] = move_agent(circleobjin,indexjes)
        agents = zeros(3,2);
        for i = 1:n
            if(indexjes(i) == 0)
                agents(i,:) = init_c_agent(i,:);
            else
                agents(i,:) = center(indexjes(i),:);
            end
        end
        delete(circleobjin);
        circleobjout = viscircles(agents,radius*ones(n,1));
        cur_agent_pos = agents;
    end
    function [center_index] = assign_agent(Cmax)
        fert_agent = find(Cmax(2,:)>0); % number of fertilizer agents
        irrig_agent = find(Cmax(1,:)>0); % number of irrigation agents
        n_fert = find(data.U_fert(:,k)>0); % indices of field being fertizlied this day
        n_irrig = find(data.U_irrig(:,k)>0); % indices of fields being iriigated this day
        if isempty(n_fert)
            n_fert = zeros(1,q_sim*length(fert_agent));
        else
            nf = q_sim*length(fert_agent)-length(n_fert);
            n_fert = [n_fert; zeros(nf,1)];
        end
        if isempty(n_irrig)
            n_irrig = zeros(q_sim*length(irrig_agent),1);
        else
            nq = q_sim*length(irrig_agent)-length(n_irrig);
            n_irrig = [n_irrig; zeros(nq,1)];
        end
        center_index(fert_agent) = n_fert((q-1)*length(fert_agent)+1:q*length(fert_agent));
        center_index(irrig_agent) = n_irrig((q-1)*length(irrig_agent)+1:q*length(irrig_agent));
    end
%% plot functions
% dont keep updating legend if variable is unchanged.
    function plotdata(save)
        if(save == 0)
            [a_temp,b_temp] = size(Current_y);
            m = min(a_temp,b_temp);
            if(m==1)
                plot(Current_x(1:k),Current_y(1:k));
            else
                if(field == (m+1))
                    leg = cell(m,1);
                    for o = 1:m
                        plot(1:k,Current_y(o,1:k));
                        if(o== 1)
                            hold on;
                        end
                        leg{o} = ['field ' num2str(o)];
                    end
                    legend(leg);
                    hold off;
                else
                    plot(1:k,Current_y(field,1:k));
                    legend(['field ' num2str(field)]);
                    
                end
                if(k~=1)
                    xlim([1 k])
                end
            end
            xlabel('Day of year')
            ylabel(varname)
        else
            [a_temp,b_temp] = size(Current_y);
            m = min(a_temp,b_temp);
            % try drawnow and animated line
            hFig = figure('visible', 'off');
            if(m==1)
                plot(Current_x(1:k),Current_y(1:k));
            else
                if(field == m+1)
                    hold on;
                    leg = cell(m,1);
                    for o = 1:m
                        plot(1:k,Current_y(o,1:k));
                        leg{o} = ['field ' num2str(o)];
                    end
                    legend(leg);
                    hold off;
                else
                    plot(1:k,Current_y(field,1:k));
                    legend(['field ' num2str(field)]);
                end
                if(k~= 1)
                    xlim([1 k])
                end
            end
            xlabel('Day of year')
            ylabel(varname)
            set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
            savename = ['figures/results_' datestr(datetime,1) '_' datestr(datetime,13)];
            savename = strrep(savename,':','_');
            savefig(hFig,savename);
            close(hFig);
        end
    end
end





