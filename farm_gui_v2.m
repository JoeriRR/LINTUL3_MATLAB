function farm_gui_v2
% Simply press the run button to start the farm GUI, an initial dataset will
% be loaded.
%
% load button: load another dataset from the folder /results, example:
% (gui1.mat).
%
% save button: saves plot to the file in /figures map as a .fig file
% current timestamp and variable name is included in savename.
%
% play button:
% simulates the controller decision with a delay of 0.5 [s].
%
% slider button: can be controlled to select day when play is not active
%
% dropdown menu fields: select fieldnumber for plot
%
% dropdown menu variable: select variable for plot

%% create figure, buttons and load simulated results
close all; clc;
data = load('results/gui.mat');
names = {'DTR', 'WN','DAVTMP','DTEFF','RAIN','U_irrig','U_fert','TSUM_r','ROOTD_r','WA_r','WC_r','WCCR_r','TEXPLO_r', 'TEVAP_r','TTRAN_r','TRUNOF_r','TIRRIG_r'...
    ,'TDRAIN_r','TRAIN_r','DVS_r', 'NNI_r','SLA_r','LAI_r','NDEMTO_r','TNSOIL_r','NUPTT_r','ANLV_r','ANST_r','ANRT_r','ANSO_r','NLOSSL_r','NLOSSR_r','WLVG_r','WST_r','WSO_r','WRT_r','TAGBM_r', ...
    'NTAC_r','LUECAL_r','CUMPAR_r','GTSUM_r','WDRT_r'};
labels = {'MJ/m^2', 'm/s','degree \circ','degree \circ','mm','mm','g/m^2','Cd','m', ...
    'mm', '-','-','mm','mm','mm','mm','mm','mm','mm','-','-','m^2/g','-', 'g/m^2','g/m^2','g/m^2', 'g/m^2','g/m^2','g/m^2','g/m^2','g/m^2','g/m^2','g/m^2','g/m^2','g/m^2'...
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
f = figure('Visible','off','Position',[100,100,1350,650]);
set(f,'MenuBar','none')
set(f,'ToolBar','none')
fieldend = 670;
figend = fieldend+500+60; buttonwidth = 100; figtop = 575;
% set ui control locations
hload = uicontrol('Style','pushbutton','String','load dataset','Position',[80,610,buttonwidth,25],'Callback',{@loadbutton_Callback});
hplay = uicontrol('Style','togglebutton','String','play','Position',[figend,20,buttonwidth,25],'Callback',{@playbutton_Callback});

% buttons for plots
hsave1 = uicontrol('Style','pushbutton','String','save','Position',[figend,figtop,buttonwidth,25],'Callback', {@savebutton1_Callback});
hsave2 = uicontrol('Style','pushbutton','String','save','Position',[figend,340-25,buttonwidth,25],'Callback', {@savebutton2_Callback});
hdata1 = uicontrol('Style','popupmenu', 'String',names,'Position',[figend+10,figtop-60,buttonwidth-20,25],'Callback',{@varpopup1_Callback});
hdata2 = uicontrol('Style','popupmenu', 'String',names,'Position',[figend+10,315-60,buttonwidth-20,25],'Callback',{@varpopup2_Callback});
fdata1 = uicontrol('Style','popupmenu', 'String',fieldnames,'Position',[figend+10,figtop-120,buttonwidth-20,25],'Callback',{@fieldpopup1_Callback});
fdata2 = uicontrol('Style','popupmenu', 'String',fieldnames,'Position',[figend+10,315-120,buttonwidth-20,25],'Callback',{@fieldpopup2_Callback});
figoptions1 = uicontrol('Style','text','String','select variable:','Position',[figend, figtop-35, buttonwidth, 25],'Fontsize',9);
figoptions2 = uicontrol('Style','text','String','select variable:','Position',[figend, 315-35, buttonwidth, 25],'Fontsize',9);
figfield1 = uicontrol('Style','text','String','select field:','Position',[figend+10, figtop-90, buttonwidth, 25],'Fontsize',9,'HorizontalAlignment', 'left');
figfield2 = uicontrol('Style','text','String','select field:','Position',[figend+10, 315-90, buttonwidth, 25],'Fontsize',9,'HorizontalAlignment', 'left');


slider_day = uicontrol('Style','slider','Position',[80,20,figend-110,25],'min',1,'max',end_date,'Callback',{@sliderbutton_Callback});
set(slider_day,'Value', 1);
dayweathertext = uicontrol('Style','text','String','DOY','Position',[80+buttonwidth, 610, 300,25], 'fontsize', 11);
% axes for schemtic view, figure1/2 and the legendbar
ha = axes('Units','pixels','Position',[80,100,520,500],'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
%ha_leg = axes('Units','pixels','Position',[figend,100,60,110],'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
ha_plot1 = axes('Units','pixels','Position',[fieldend+30,360,500,240],'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
ha_plot2 = axes('Units','pixels','Position',[fieldend+30,100,500,240],'box','off','xtick',[],'ytick',[],'ztick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);
% staticlegStr = ["Wet"', "Optimal growth", "Normal growth", "No growth"]; staticlegStr = pad(staticlegStr);
% staticleg1 = uicontrol('Style','text','String',strcat(staticlegStr(1)),'Position',[figend+70, 179, buttonwidth, 25],'Fontsize',10,  'HorizontalAlignment', 'left');
% staticleg2 = uicontrol('Style','text','String',strcat(staticlegStr(2)),'Position',[figend+70, 152, buttonwidth, 25],'Fontsize',10,'HorizontalAlignment', 'left');
% staticleg3 = uicontrol('Style','text','String',strcat(staticlegStr(3)),'Position',[figend+70, 124, buttonwidth, 25],'Fontsize',10,'HorizontalAlignment', 'left');
% staticleg4 = uicontrol('Style','text','String',strcat(staticlegStr(4)),'Position',[figend+70, 98, buttonwidth, 25],'Fontsize',10,'HorizontalAlignment', 'left');
f.Visible = 'on';

%% initialise UI and make all buttons resizable
f.Units = 'normalized';
%ha_leg.Units = 'normalized';
ha.Units = 'normalized';
ha_plot1.Units = 'Normalized';
ha_plot2.Units = 'Normalized';
hplay.Units = 'normalized';
dayweathertext.Units = 'Normalized';
hload.Units = 'normalized';
fdata1.Units = 'normalized';
fdata2.Units = 'normalized';

figoptions1.Units = 'normalized';
figoptions2.Units = 'Normalized';
figfield1.Units = 'normalized';
figfield2.Units = 'normalized';
slider_day.Units = 'normalized';
hdata1.Units = 'normalized';
hdata2.Units = 'normalized';

hsave1.Units = 'normalized';
hsave2.Units = 'normalized';
% staticleg1.Units = 'normalized';
% staticleg2.Units = 'normalized';
% staticleg3.Units = 'normalized';
% staticleg4.Units = 'normalized';
range = [-agentlength farmlength 0 farmlength];
axis(range);
global Current_x Current_y1 Current_y2 label varname field k q data_changed ax_in;

%% the GUI loop
while(true)
    if(isvalid(f))
        switch SwitchCase
            case 0 % init
                % draw_legend(figend,buttonwidth,ha_leg)
                [init_c_agent,radius,rect_set,center,~,circles] = init_after_data_loaded; %initiase uifigire and variables
                cur_agent_pos = init_c_agent;                
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
                % check if data changed
                if(data_changed == true) % only plot when data is changed
                    ax_in = [ha_plot1 ha_plot2];
                    plotdata(0,ax_in); % plot data
                end
                data_changed = false; % set data changed false, changed in buttonpresses to true
                update_field_color(rect_set);% redraw initial field at current day
                [indecis,resources,fert_agent,irrig_agent] = assign_agent(data.Cmax);
                [circles] = move_agent(circles,indecis,fert_agent,irrig_agent,init_c_agent,radius);
                res = data.Cmax-resources; % if no refil required substract from previous in subtime
                set_Day_Weather(data); % always plot dayweather data
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
    function varpopup1_Callback(hdata,source,eventdata)
        val = hdata.Value; % reads selected variable used for plot
        Str = hdata.String;
        varname{1} = Str{val};
        label{1} = labels{val};
        Current_y1 = data.(varname{1}); % select current variable from loaded structure
        data_changed =true;
    end
    function varpopup2_Callback(hdata,source,eventdata)
        val = hdata.Value; % reads selected variable used for plot
        Str = hdata.String;
        varname{2} = Str{val};
        label{2} = labels{val};
        Current_y2 = data.(varname{2}); % select current variable from loaded structure
        data_changed =true;
    end
    function fieldpopup1_Callback(fdata,source,eventdata)
        val = fdata.Value; % reads fielddata to plot
        Str = fdata.String;
        if(strcmp(Str{val},'all')) % used for all plot
            field(1) = m+1;
        else
            field(1) = str2double(Str{val});
        end
        data_changed = true;
    end
    function fieldpopup2_Callback(fdata,source,eventdata)
        val = fdata.Value; % reads fielddata to plot
        Str = fdata.String;
        if(strcmp(Str{val},'all')) % used for all plot
            field(2) = m+1;
        else
            field(2) = str2double(Str{val});
        end
        data_changed = true;
    end
    function sliderbutton_Callback(hplay,source,eventdata)
        data_changed = true; % if plot is active and slider is used, data will have to be updated
    end
    function savebutton1_Callback(source,eventdata)
        % it plots the same data but this time not in the ui window to
        % save the plot only and not the whole ui figure
        % create invisible figure, save it to a visibile .fig file and
        % remove the figure
        plotdata(1,ax_in);
    end
    function savebutton2_Callback(source,eventdata)
        % it plots the same data but this time not in the ui window to
        % save the plot only and not the whole ui figure
        % create invisible figure, save it to a visibile .fig file and
        % remove the figure
        plotdata(2,ax_in);
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
%{
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
%}
%{
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
%}
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
                text(center(index,1),center(index,2),num2str(index),'Color',[1 1 1],'HorizontalAlignment','center');
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
    function [circleobjout] = move_agent(circleobjin,indexjes,fert_agent,irrig_agent,init_c_agent,radius)
        agents = zeros(3,2);
        for i = 1:n
            if(indexjes(i) == 0 || indexjes(i) > m)
                % agents(i,:) = init_c_agent(i,:);
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
        if(isfield(data,'A') == 0) % change this since it does not work for old versions
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
        else
            actual_res = zeros(size(Cmax));
            fert_agent = find(Cmax(2,:)>0);
            irrig_agent = find(Cmax(1,:)>0);
            for i = 1:n
               center_index(i) = find(data.A{k}(i,:,q)>0);                
            end
        end
        
    end
    function update_field_color(rect_set)
        % add these variables to gui.mat create colorbar
        WCWP = 0.2;
        WCFC = 0.4;
        rgbvalues(1,:) = [102 178 255];
        rgbvalues(2,:) = [51 153 255];
        rgbvalues(3,:) = [0 128 255];
        rgbvalues(4,:) = [0 0 255];
        for i = 1:length(rect_set)
            if(data.WC_r(i,k)<WCWP)
                rect_set(i).FaceColor = [1 0 0];
            end
            if(data.WC_r(i,k)>=WCWP && data.WC_r(i,k)<data.WCCR_r(i,k))
                rect_set(i).FaceColor = [1 1 0];
            end
            if(data.WC_r(i,k)>= data.WCCR_r(i,k) && data.WC_r(i,k)<=WCFC)
                percentage = (1-(WCFC-data.WC_r(i,k))/(WCFC-data.WCCR_r(i,k)))*4;
                for z = 1:length(rgbvalues)
                    if(percentage<(z))
                        percentage = z;
                        break;
                    end
                end
                rect_set(i).FaceColor = rgbvalues(percentage,:)/255;
            end
            if(data.WC_r(i,k)>WCFC)
                rect_set(i).FaceColor = [0 0 0.5];
            end
        end
        
    end
%{
function draw_legend(figend,buttonwidth,ha_leg)
        axes(ha_leg)
        axis([0,10,0,10])
        rectangle('Position',[0,0,10,2.5],'FaceColor',[1 1 0]);
        rectangle('Position',[0,2.5,10,5],'FaceColor',[0.5843 0.8157 0.9882]);
        rectangle('Position',[0,5,10,7.5],'FaceColor',[0 0 1]);
        rectangle('Position',[0,7.5,10,10],'FaceColor',[0 0 0.5]);
        axes(ha)
    end
%}
%% plot functions
    function plotdata(save,ax_in)
        if(save == 0)
            for figcount = 1:length(ax_in)
                if(figcount == 1)
                    Current_y = Current_y1;
                else
                    Current_y = Current_y2;
                end
                axes(ax_in(figcount));
                [a_temp,b_temp] = size(Current_y);
                m = min(a_temp,b_temp);
                tempvar = strrep(varname{figcount},'_r', '');
                tempvar = strrep(tempvar,'_', '');
                if(m==1)
                    plot(Current_x(1:k),Current_y(1:k));
                    title([tempvar, ' for field ' num2str(field(figcount))]);
                else
                    if(field(figcount) == (m+1))
                        for o = 1:m
                            plot(Current_x(1:k),Current_y(o,1:k));
                            if(o==1)
                                hold on;
                            end
                        end
                        title([tempvar,' all fields']);
                        hold off;
                    else
                        plot(Current_x(1:k),Current_y(field(figcount),1:k));
                        title([tempvar, ' for field ' num2str(field(figcount))]);
                    end
                end
                if(k~=1)
                    xlim([Current_x(1) Current_x(k)])
                end
                if(figcount == 2)
                    xlabel('Day of year')
                else
                    set(ax_in(figcount),'Xticklabel',[]);
                end
                ylabel(label{figcount})
            end
        else
            if(save == 1)
                Current_y = Current_y1;
            else
                Current_y = Current_y2;
            end
            hFig = figure('visible', 'off');
            [a_temp,b_temp] = size(Current_y);
            m = min(a_temp,b_temp);
            tempvar = strrep(varname{save},'_r', '');
            tempvar = strrep(tempvar,'_', '');
            if(m==1)
                plot(Current_x(1:k),Current_y(1:k));
                title([tempvar, ' for field ' num2str(field(save))]);
            else
                if(field(save) == (m+1))
                    for o = 1:m
                        plot(Current_x(1:k),Current_y(o,1:k));
                        if(o==1)
                            hold on;
                        end
                    end
                    title([tempvar,' all fields']);
                    hold off;
                else
                    plot(Current_x(1:k),Current_y(field(save),1:k));
                    title([tempvar, ' for field ' num2str(field(save))]);
                end
            end
            if(k~=1)
                xlim([Current_x(1) Current_x(k)])
            end
            xlabel('Day of year')
            ylabel(label{save})
            set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
            savename = ['figures/' varname{save} '_' datestr(datetime,1) '_' datestr(datetime,13)];
            savename = strrep(savename,':','_');
            savefig(hFig,savename);
            close(hFig);
        end
        
        axes(ha);
    end
%% initialise functions
    function [init_c_agent,radius,rect_set,center,center_agent,circles] = init_after_data_loaded
        m = size(data.LAI_r,1); % get number of fields m = 10
        [r,n] = size(data.Cmax); % get number of resources
        for bb=1:m
            fieldnames{bb} = num2str(bb); % write number of fields to cell array string for button control
        end
        fieldnames{bb+1} = 'all'; % posibility to plot all fields at once
        [end_date] = length(data.U_irrig); q_sim = data.q; % simulation stop date and subtime units
        % reset data dependent ui controls
        fdata1.Value = 1;
        fdata2.Value = 1;
        fdata1.String = fieldnames;
        fdata2.String = fieldnames;
        slider_day.Max = end_date;
        slider_day.Value = 1;
        
        % set variable to initial state
        q = 1; % subtime unit
        k = 1;
        % draw field and agents at initial stage
        axes(ha);
        cla;
        [center, fieldwidth,fieldheigth,rect_set]= draw_Field(m,range,5/7);
        update_field_color(rect_set);
        radius = 1/4*min(fieldwidth,fieldheigth);
        [circles,init_c_agent] = init_Agent(n,range,radius);
        center_agent = init_c_agent; draw_Border(range);
        
        % print day and weather status
        set_Day_Weather(data);
        
        % initialise plot variables
        Current_x = data.DOY; Current_y1 = data.RAIN; Current_y2 = data.LAI_r; varname{1} = names{5}; varname{2} = names{23}; label{1} = labels{5}; label{2} = labels{23};
        field = ones(2,1); data_changed = true;
        
        % go to loop
        SwitchCase = 1;
    end
end





