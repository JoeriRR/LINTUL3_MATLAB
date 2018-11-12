function farm_gui
% to do % dtr wn rain
%% create figure, buttons and load simulated results
close all; clc;
data = load('results/gui.mat');
names = {'DTR', 'WN','DAVTMP','DTEFF','RAIN','U_irrig','U_fert','TSUM_r','ROOTD_r','WA_r','WC_r','WCCR_r','TEXPLO_r', 'TEVAP_r','TTRAN_r','TRUNOF_r','TIRRIG_r'...
    ,'TDRAIN_r','TRAIN_r','DVS_r', 'NNI_r','SLA_r','LAI_r','NDEMTO_r','TNSOIL_r','NUPTT_r','ANLV_r','ANST_r','ANRT_r','ANSO_r','NLOSSL_r','NLOSSR_r','WLVG_r','WST_r','WSO_r','WRT_r','TAGBM_r', ...
    'NTAC_r','LUECAL_r','CUMPAR_r','GTSUM_r','WDRT_r'};
labels = {'MJ/m^2', '?','degree \circ','degree \circ','mm','mm','gN/m^2','Cd','m', ...
    'mm', '-','-','mm','mm','mm','mm','mm','mm','mm','-','-','m^2/g','-', 'gN/m^2','gN/m^2','gN/m^2', 'g/m^2','g/m^2','g/m^2','g/m^2','gN/m^2','gN/m^2','g/m^2','g/m^2','g/m^2'...
    ,'g/m^2','g/m^2','g/m^2','MJ/m^2','g/m^2','g/m^2','g/m^2'};
% initiasition of some variables
farmlength = 10;
agentlength = 2;
SwitchCase = 0;
% asssign varibles ie qmax,m,n,radius etc. create function for this
m = size(data.LAI_r,1); 
[r,n] = size(data.Cmax);
fieldnames = cell(1,m);
for b=1:m
    fieldnames{b} = num2str(b);
end
fieldnames{b+1} = 'all';
[end_date] = length(data.U_irrig); q_sim = data.q; % need to add U and delta matrix to assign the agents!
%% buttons and text
%position [left bottom width heigth]!
f = figure('Visible','off','Position',[100,200,900,600]);
set(f,'MenuBar','none')
set(f,'ToolBar','none')
figend = 670; buttonwidth = 100; figtop = 525;
% set ui control locations
hsave = uicontrol('Style','pushbutton','String','save','Position',[figend+110,figtop-190,buttonwidth,25],'Callback', {@savebutton_Callback});
hload = uicontrol('Style','pushbutton','String','load','Position',[figend,figtop-190,buttonwidth,25],'Callback',{@loadbutton_Callback});
hplay = uicontrol('Style','togglebutton','String','play','Position',[figend,figtop,buttonwidth,25],'Callback',{@playbutton_Callback});
hplot = uicontrol('Style','togglebutton','String','plot','Position',[figend+110,figtop,buttonwidth,25],'Callback',{@plotbutton_Callback});
hdata = uicontrol('Style','popupmenu', 'String',names,'Position',[figend+10,figtop-60,buttonwidth-20,25],'Callback',{@varpopup_Callback});
fdata = uicontrol('Style','popupmenu', 'String',fieldnames,'Position',[figend+120,figtop-60,buttonwidth-20,25],'Callback',{@fieldpopup_Callback});
figoptions = uicontrol('Style','text','String','select variable:','Position',[figend, figtop-35, buttonwidth, 25],'Fontsize',9);
figfield = uicontrol('Style','text','String','select field:','Position',[figend+110, figtop-35, buttonwidth, 25],'Fontsize',9);
slider_day = uicontrol('Style','slider','Position',[80,20,550,25],'min',1,'max',end_date,'Callback',{@sliderbutton_Callback});
set(slider_day,'Value', 1);
StaticField = uicontrol('Style','text', 'String','field status:','Position', [figend, figtop-80, 100, 15],'Fontsize',9);
StaticField2 = uicontrol('Style','text', 'String','agent status:','Position', [figend+110, figtop-80, 100, 15],'Fontsize',9);
dayweathertext = uicontrol('Style','text','String','DOY','Position',[80, 550, 300,25], 'fontsize', 11);
fieldtext = uicontrol('Style','text','String','-','Position',[figend, figtop-140, 100, 50]);
agenttext = uicontrol('Style','text','String','-','Position',[figend+110, figtop-140, 100, 50]);
ha = axes('Units','pixels','Position',[80,100,550,450],'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
ha_leg = axes('Units','pixels','Position',[figend,100,60,110],'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);

staticlegStr = ["Wet"', "Optimal growth", "Normal growth", "No growth"]; staticlegStr = pad(staticlegStr);
staticleg1 = uicontrol('Style','text','String',strcat(staticlegStr(1)),'Position',[figend+70, 179, buttonwidth, 25],'Fontsize',10,  'HorizontalAlignment', 'left');
staticleg2 = uicontrol('Style','text','String',strcat(staticlegStr(2)),'Position',[figend+70, 152, buttonwidth, 25],'Fontsize',10,'HorizontalAlignment', 'left');
staticleg3 = uicontrol('Style','text','String',strcat(staticlegStr(3)),'Position',[figend+70, 124, buttonwidth, 25],'Fontsize',10,'HorizontalAlignment', 'left');
staticleg4 = uicontrol('Style','text','String',strcat(staticlegStr(4)),'Position',[figend+70, 98, buttonwidth, 25],'Fontsize',10,'HorizontalAlignment', 'left');
f.Visible = 'on';

%% initialise UI and make all buttons resizable
f.Units = 'normalized';
ha_leg.Units = 'normalized';
ha.Units = 'normalized';
hplay.Units = 'normalized';
fieldtext.Units = 'normalized';
agenttext.Units = 'normalized';
StaticField.Units = 'normalized';
StaticField2.Units = 'normalized';
dayweathertext.Units = 'Normalized';
hload.Units = 'normalized';
hplot.Units = 'normalized';
fdata.Units = 'normalized';
figoptions.Units = 'normalized';
figfield.Units = 'normalized';
slider_day.Units = 'normalized';
hdata.Units = 'normalized';
hsave.Units = 'normalized';
staticleg1.Units = 'normalized';
staticleg2.Units = 'normalized';
staticleg3.Units = 'normalized';
staticleg4.Units = 'normalized';
range = [-agentlength farmlength 0 farmlength];
axis(range); label = labels{23}; axes(ha); current_field = []; current_agent = [];
%% the GUI loop
while(true)
    if(isvalid(f))
        switch SwitchCase
            case 0 % init
                draw_legend(figend,buttonwidth,ha_leg)
                [Current_x,Current_y,field,init_c_agent,radius,varname] = init_after_data_loaded; %initiase uifigire and variables
            case 1 % main loop
                if(get(hplay,'Value')) % start incrementing day if play button is active                   
                    if (k >= end_date)
                        k = end_date; % stop when endate reached
                    else
                        if(mod(q,q_sim) == 0) % update day when q subtime units have passed
                            q = 0;
                            k = k+1;
                            data_changed = true;
                        end
                        q = q+1;
                    end
                    set(slider_day,'value',k); % if not active, slider control can be used to select day
                else
                    k = round(get(slider_day,'value')); % round since day is a int
                end
                if(get(hplot,'Value')) % plot variables
                    if(last_entry == true)
                        iptremovecallback(f,'WindowButtonMotionFcn', id1); % revove mouse hover function to see field/agent status
                    end
                    % check if data changed
                    if(data_changed == true) % only plot when data is changed
                        cla;
                        plotdata(0); % plot data
                    end
                    data_changed = false; % set data changed false, changed in buttonpresses to true
                    last_entry = false; 
                else
                    % only do once
                    if(last_entry == false) % add button listener
                        id1 = iptaddcallback(f,'WindowButtonMotionFcn',@selecter);
                        axis off; % remove axis generated by figure
                        cla;
                        legend('hide') % hide legend if still active 
                        [center, fieldwidth,fieldheigth,rect_set]= draw_Field(m,range,5/7);
                        update_field_color(rect_set);% redraw initial field at current day
                        [indecis,resources,fert_agent,irrig_agent] = assign_agent(data.Cmax);
                        [circles,center_agent] = move_agent(circles,indecis,fert_agent,irrig_agent);
                        draw_Border(range);
                        axis(range); % not sure if required
                    else
                        [indecis,resources,fert_agent,irrig_agent] = assign_agent(data.Cmax); % reads the required resources                        
                        [circles,center_agent] = move_agent(circles,indecis,fert_agent,irrig_agent); % move agents to current assigned field in subtime unit q
                        update_field_color(rect_set);% redraw initial field at current day

                    end
                    res = data.Cmax-resources; % if no refil required substract from previous in subtime
                    last_entry = true;
                end
                set_Day_Weather(data); % always plot dayweather data
                update_selected;
                pause(0.5); % pause the gui for some time, it will be laggy when no pause is used
        end
    else
        return; % close while loop when exit button is pressed
    end
end
%% button functions
    function playbutton_Callback(source,eventdata)
        if(get(hplay,'Value')) % change names depending on the active state
            set(hplay,'String', 'pause')
        else
            set(hplay,'String', 'play')
        end
    end
    function plotbutton_Callback(source,eventdata)
        if(get(hplot,'Value') == 1) % change names depending on the active state
            set(hplot,'String', 'vis')
        else
            set(hplot,'String','plot')
        end
        data_changed = true;
    end
    function varpopup_Callback(hdata,source,eventdata) 
        val = hdata.Value; % reads selected variable used for plot       
        Str = hdata.String; 
        varname = Str{val};
        label = labels{val};
        Current_y = data.(varname); % select current variable from loaded structure
        data_changed =true;
    end
    function fieldpopup_Callback(fdata,source,eventdata)
        val = fdata.Value; % reads fielddata to plot 
        Str = fdata.String;
        if(strcmp(Str{val},'all')) % used for all plot 
            field = m+1;
        else
            field = str2double(Str{val});
        end
        data_changed = true;
    end
    function sliderbutton_Callback(hplay,source,eventdata)
        data_changed = true; % if plot is active and slider is used, data will have to be updated
    end
    function savebutton_Callback(source,eventdata)
        if(get(hplot,'value'))
            % it plots the same data but this time not in the ui window to
            % save the plot only and not the whole ui figure 
            % create invisible figure, save it to a visibile .fig file and
            % remove the figure
            plotdata(1);
        end
    end
    function loadbutton_Callback(source,eventdata)
        path = 'results/'; path = fullfile(path,'*.mat');
        [Filename,Pathname] = uigetfile(path,'Select mat file');
        if Filename == 0
           return; 
        end
        fullFileName = fullfile(Pathname,Filename);
        clearvars data fieldnames;
        data = load(fullFileName);
        k = 1;
        SwitchCase = 0;
    end
%% textbox functions
    function set_Day_Weather(data)
        DOY = data.DOY;  % print of q, day and current rainfall.
        RAIN = data.RAIN;
        set(dayweathertext,'String',['day ' num2str(DOY(k)) ' subtime unit q ' num2str(q) ' rainfall ' num2str(RAIN(k)) '[mm]']);
    end
    function [current_pos] = selecter(source,eventdata)
        % ad previous selecter
        current_pos = get(gca,'CurrentPoint'); % check whether mouse in one of the fields
        for i = 1:size(center,1)
            if(current_pos(1,1)>=(center(i,1)-fieldwidth/2) && current_pos(1,1)<= (center(i,1)+fieldwidth/2) && current_pos(1,2)>=(center(i,2)-fieldheigth/2) && current_pos(1,2)<= (center(i,2)+fieldheigth/2))
                set(fieldtext,'String',['Field ',num2str(i), newline, ['LAI ' num2str(round(data.LAI_r(i,k),1)) ' [-]'], newline, ['WSO ' num2str(round(data.WSO_r(i,k),1)) ' [g/m^2]']]);
                current_field = i;
                break;  % leave loop since mouse can only be in at most on field         
            end
        end
        for i = 1:size(center_agent,1) % check whether mouse in inside on of the agents
            r = sqrt((current_pos(1,1)-center_agent(i,1))^2+(current_pos(1,2)-center_agent(i,2))^2);
            if(r<=radius)
                set(agenttext,'String',['Agent ',num2str(i), newline, ['water ' num2str(res(1,i)) ' [mm]'], newline, ['nitrogen ' num2str(res(2,i))] ' [gN/m^2]']);
                current_agent = i;
                break; % leave loop since mouse can only be in at most on agent         
            end
        end
    end
    function update_selected
        if(isempty(current_field)~=1)
            set(fieldtext,'String',['Field ',num2str(current_field), newline, ['LAI ' num2str(round(data.LAI_r(current_field,k),1)) ' [-]'], newline, ['WSO ' num2str(round(data.WSO_r(current_field,k),1)) ' [g/m^2]']]);
        else
            set(fieldtext,'String','');
        end
        if(isempty(current_agent)~=1)
            set(agenttext,'String',['Agent ',num2str(current_agent), newline, ['water ' num2str(res(1,current_agent)) ' [mm]'], newline, ['nitrogen ' num2str(res(2,current_agent))] ' [gN/m^2]']);
        else
            set(agenttext,'String','');
        end
    end
%% draw functions
    function [center,fieldwidth,fieldheigth,rect_set] = draw_Field(m,range,fillfact)
        % draw fields as rectangles
        %figure(f);
        rows = round(sqrt(m));
        if (mod(m,rows) == 0)
            columns = m/rows;
        else
            if(rows^2<m)
                columns = rows+1;
            else
                columns = rows;
            end
        end
        fieldwidth = abs(0-range(2))*fillfact/columns;
        openwidth = (abs(0-range(2))-fieldwidth*columns)/(columns+1);
        fieldheigth = abs(range(3)-range(4))*fillfact/rows;
        openheigth = (abs(range(3)-range(4))-fieldheigth*rows)/(rows+1);
        index = 1; center = zeros(m,2);
        for i = 1:rows
            for j = 1:columns
                rect_set(index) = rectangle('Position',[openwidth*j+(j-1)*fieldwidth openheigth*i+(i-1)*fieldheigth fieldwidth fieldheigth]);
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
    function [circleobj,center] = init_Agent(n,range,radius)
        index = 1;
        gap = (abs(range(3)-range(4))-radius*n*2)/(n+1);
        center = zeros(n,2);
        for i = 1:n
            center(index,1:2) = [range(1)/2, gap*i+(radius*(2*i-1))];
            index = index+1;
        end
        circleobj = viscircles(center,radius*ones(n,1));
    end
    function [circleobjout,cur_agent_pos] = move_agent(circleobjin,indexjes,fert_agent,irrig_agent)
        agents = zeros(3,2);
        for i = 1:n
            if(indexjes(i) == 0)
                agents(i,:) = init_c_agent(i,:);
            else
                agents(i,:) = center(indexjes(i),:);
            end
        end
        delete(circleobjin);
        circleobjout(1) = viscircles(agents(irrig_agent,:),radius*ones(length(irrig_agent),1),'Color','green');
        circleobjout(2) = viscircles(agents(fert_agent,:),radius*ones(length(fert_agent),1),'Color',[165/255 42/255 42/255]);
        cur_agent_pos = agents;
    end
    function [center_index,actual_res,fert_agent,irrig_agent] = assign_agent(Cmax)
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
        actual_res = zeros(size(Cmax));
        for i = 1:n
            if(center_index(i) ~= 0)
                if(i<=length(irrig_agent))
                    actual_res(1,i) = data.U_irrig(center_index(i),k);
                else
                    actual_res(2,i) = data.U_fert(center_index(i),k);
                end
            end
        end
    end
    function update_field_color(rect_set)
       % add these variables to gui.mat
       WCWP = 0.2;
       WCFC = 0.4;       
       for i = 1:length(rect_set)
           if(data.WC_r(i,k)<WCWP)
               rect_set(i).FaceColor = [1 1 0];
           end
           if(data.WC_r(i,k)>=WCWP && data.WC_r(i,k)<data.WCCR_r(i,k))
               rect_set(i).FaceColor = [0.5843 0.8157 0.9882];
           end
           if(data.WC_r(i,k)>= data.WCCR_r(i,k) && data.WC_r(i,k)<=WCFC)
               rect_set(i).FaceColor = [0 0 1];
           else
               rect_set(i).FaceColor = [0 0 0.5];
           end           
       end       
       
    end
    function draw_legend(figend,buttonwidth,ha_leg)
        axes(ha_leg)
        axis([0,10,0,10])
        rectangle('Position',[0,0,10,2.5],'FaceColor',[1 1 0]);
        rectangle('Position',[0,2.5,10,5],'FaceColor',[0.5843 0.8157 0.9882]);
        rectangle('Position',[0,5,10,7.5],'FaceColor',[0 0 1]);
        rectangle('Position',[0,7.5,10,10],'FaceColor',[0 0 0.5]);  
        axes(ha)
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
                        plot(Current_x(1:k),Current_y(o,1:k));
                        if(o== 1)
                            hold on;
                        end
                        leg{o} = ['field ' num2str(o)];
                    end
                    legend(leg);
                    hold off;
                else
                    plot(Current_x(1:k),Current_y(field,1:k));
                    legend(['field ' num2str(field)]);
                    
                end
                if(k~=1)
                    xlim([Current_x(1) Current_x(k)])
                end
            end
            xlabel('Day of year')
            ylabel(label)
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
                        plot(Current_x(1:k),Current_y(o,1:k));
                        leg{o} = ['field ' num2str(o)];
                    end
                    legend(leg);
                    hold off;
                else
                    plot(Current_x(1:k),Current_y(field,1:k));
                    legend(['field ' num2str(field)]);
                end
                if(k~= 1)
                    xlim([Current_x(1) Current_x(k)])
                end
            end
            xlabel('Day of year')
            ylabel(label)
            set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
            savename = ['figures/' varname '_' datestr(datetime,1) '_' datestr(datetime,13)];
            savename = strrep(savename,':','_');
            savefig(hFig,savename);
            close(hFig);
        end
    end
%% initialise functions
    function [Current_x,Current_y,field,init_c_agent,radius,varname,rect_set] = init_after_data_loaded
        m = size(data.LAI_r,1); % get number of fields m = 10
        [r,n] = size(data.Cmax); % get number of resources
        for bb=1:m
            fieldnames{bb} = num2str(bb); % write number of fields to cell array string for button control
        end
        fieldnames{bb+1} = 'all'; % posibility to plot all fields at once
        [end_date] = length(data.U_irrig); q_sim = data.q; % simulation stop date and subtime units
        % reset data dependent ui controls
        fdata.String = fieldnames;
        slider_day.Max = end_date;
        slider_day.Value = 1;
        
        % set variable to initial state       
        q = 1; % subtime unit
        k = 1;
        % draw field and agents at initial stage
        [center, fieldwidth,fieldheigth,rect_set]= draw_Field(m,range,5/7);
        update_field_color(rect_set);
        radius = 1/4*min(fieldwidth,fieldheigth);
        [circles,init_c_agent] = init_Agent(n,range,radius);
        center_agent = init_c_agent; draw_Border(range);
        
        % print day and weather status
        set_Day_Weather(data);
        id1 = iptaddcallback(f,'WindowButtonMotionFcn',@selecter);
        
        % initialise plot variables
        Current_x = data.DOY; Current_y = data.LAI_r; varname ='LAI';
        field = 1; last_entry = false; data_changed = true;
        
        % go to loop  
        current_agent = []; current_field = [];
        SwitchCase = 1;        
    end
end





