function ax=ggplotAxes2D(varargin)
%
% @author:slandarer
% 
%  Parameter description ?
% -----------------------------------------------------
% AxesTheme   |  Coordinate area style        | 'gray'/'economist'/'wsj'/'own1'
% ColorOrder  |  Drawing object color sequence    | 'default'/'none'/'npg'/'lancet'/'starterk'
%?????????????????  'Set1'/'Set2'/'Set3'/'Dark2'/'own1'
% LegendStyle |  Legend style            | 'ggplot'/'own1'
% EdgeStyle   |  Outline style            | 'none'/'gray'/'white'/'ori'

% ax.Legend.UserData.NewBkg    Legend new background 
% ax.Legend.UserData.NewTitle  Legend new title 
% ax.UserData.NewYTick(i)     Y Axis new label 

%  Get the coordinate area to process =====================================================
if strcmp(get(varargin{1},'type'),'axes' )
    ax=varargin{1};
else
    ax=gca;
end
hold(ax,'on')



% default==================================================================
theme.AxesTheme='gray';
theme.ColorOrder='default';
theme.LegendStyle='ggplot';
theme.EdgeStyle='none';

% Extract useful information from variable length variables ==============================================
for i=1:length(varargin)
    tempVar=varargin{i};
    if strcmp(tempVar,'AxesTheme')||strcmp(tempVar,'axesTheme')||strcmp(tempVar,'axestheme')
        theme.AxesTheme=varargin{i+1};
    end
    if strcmp(tempVar,'ColorOrder')||strcmp(tempVar,'colorOrder')||strcmp(tempVar,'colororder')
        theme.ColorOrder=varargin{i+1};
    end
    if strcmp(tempVar,'LegendStyle')||strcmp(tempVar,'legendStyle')||strcmp(tempVar,'legendstyle')
        theme.LegendStyle=varargin{i+1};
    end
    if strcmp(tempVar,'EdgeStyle')||strcmp(tempVar,'edgeStyle')||strcmp(tempVar,'edgestyle')
        theme.EdgeStyle=varargin{i+1};
    end
end


%  The color scheme 
switch theme.ColorOrder
    case 'none'
    case 'default'
        ax.ColorOrder=[0.9900    0.4500    0.4500
    0.8500    0.5600         0
    0.6400    0.6500         0
    0.2200    0.7100         0
         0    0.7500    0.4900
         0    0.7500    0.7700
         0    0.6900    0.9600
    0.5800    0.5600    1.0000
    0.9100    0.4200    0.9500
    1.0000    0.3800    0.7400];   
    case 'npg'
        ax.ColorOrder=[0.9000    0.2900    0.2100
    0.3000    0.7300    0.8400
         0    0.6300    0.5300
    0.2400    0.3300    0.5300
    0.9500    0.6100    0.5000
    0.5200    0.5700    0.7100
    0.5700    0.8200    0.7600
    0.8600         0         0
    0.4900    0.3800    0.2800
    0.6900    0.6100    0.5200];
    case 'lancet'
        ax.ColorOrder=[     0    0.2700    0.5500
    0.9300         0         0
    0.2600    0.7100    0.2500
         0    0.6000    0.7100
    0.5700    0.3700    0.6200
    0.9900    0.6900    0.5700
    0.6800         0    0.1600
    0.6800    0.7100    0.7100
    0.1100    0.1000    0.1000];
    case 'starterk'
        ax.ColorOrder=[0.8000    0.0500         0
    0.3600    0.5300    0.8500
    0.5200    0.7400         0
    1.0000    0.8000         0
    0.4900    0.5300    0.5600
         0    0.7100    0.8900
         0    0.6900    0.4000];
    case 'Set1'
        ax.ColorOrder=[0.8900    0.1000    0.1100
    0.2200    0.4900    0.7200
    0.3000    0.6900    0.2900
    0.6000    0.3100    0.6400
    1.0000    0.5000         0
    1.0000    1.0000    0.2000
    0.6500    0.3400    0.1600
    0.9700    0.5100    0.7500
    0.6000    0.6000    0.6000];
    case 'Set2'
        ax.ColorOrder=[0.4000    0.7600    0.6500
    0.9900    0.5500    0.3800
    0.5500    0.6300    0.8000
    0.9100    0.5400    0.7600
    0.6500    0.8500    0.3300
    1.0000    0.8500    0.1800
    0.9000    0.7700    0.5800
    0.7000    0.7000    0.7000];
    case 'Set3'
        ax.ColorOrder=[0.5500    0.8300    0.7800
    1.0000    1.0000    0.7000
    0.7500    0.7300    0.8500
    0.9800    0.5000    0.4500
    0.5000    0.6900    0.8300
    0.9900    0.7100    0.3800
    0.7000    0.8700    0.4100
    0.9900    0.8000    0.9000
    0.8500    0.8500    0.8500
    0.7400    0.5000    0.7400
    0.8000    0.9200    0.7700
    0.8300    0.8300    0.8300];
    case 'Dark2'
        ax.ColorOrder=[0.1100    0.6200    0.4700
    0.8500    0.3700    0.0100
    0.4600    0.4400    0.7000
    0.9100    0.1600    0.5400
    0.4000    0.6500    0.1200
    0.9000    0.6700    0.0100
    0.6500    0.4600    0.1100
    0.4000    0.4000    0.4000];
    case 'own1'
        ax.ColorOrder=[0.8500    0.7100    0.8000
    0.3700    0.4400    0.6600
    0.7500    0.6900    0.8300
    0.3700    0.2200    0.5200
    0.8400    0.2500    0.5500
    0.7200    0.5200    0.5200
    0.6100    0.3800    0.6000
    0.0400    0.1400    0.2800
    1.0000    0.5800    0.3500
    0.9500    0.8900    0.7500];
end

%  part plot scatter modification 
if false
childrenNum=length(ax.Children);
for i=1:childrenNum
    switch theme.EdgeStyle
        case 'none'
            EdgeColor=[];
        case 'gray'
            EdgeColor=[0.3 0.3 0.3];
        case 'white'
            EdgeColor=[0.96 0.96 0.96];
        case 'ori'
            EdgeColor=ax.ColorOrder(mod(i-1,size(ax.ColorOrder,1))+1,:).*0.5;
    end
    switch get(ax.Children(i),'type')
        case 'line'
            ax.Children(i).LineWidth=1.8;
            if ~isempty(EdgeColor)
                ax.Children(i).LineWidth=1.5;
                ax.Children(i).MarkerEdgeColor=EdgeColor;
                ax.Children(i).MarkerSize=8;
                ax.Children(i).MarkerFaceColor=ax.ColorOrder(mod(i-1,size(ax.ColorOrder,1))+1,:);
            end
        case 'scatter'
            if ~isempty(EdgeColor)
                ax.Children(i).MarkerFaceColor=ax.ColorOrder(mod(i-1,size(ax.ColorOrder,1))+1,:);
                ax.Children(i).LineWidth=1.5;
                ax.Children(i).MarkerEdgeColor=EdgeColor;
                ax.Children(i).SizeData=60;
            end
    end
end
end


if true
childrenNum=length(ax.Children);
lineSet=[];
scatterSet=[];
n=1;
for i=childrenNum:-1:1
    if strcmp(get(ax.Children(i),'type'),'scatter')
        scatterSet{n}=ax.Children(i);
        n=n+1;
    end
end
n=1;
for i=childrenNum:-1:1
    if strcmp(get(ax.Children(i),'type'),'line')
        lineSet{n}=ax.Children(i);
        n=n+1;
    end
end
for i=1:length(scatterSet)
    switch theme.EdgeStyle
        case 'none'
            EdgeColor=[];
        case 'gray'
            EdgeColor=[0.3 0.3 0.3];
        case 'white'
            EdgeColor=[0.96 0.96 0.96];
        case 'ori'
            EdgeColor=ax.ColorOrder(mod(i-1,size(ax.ColorOrder,1))+1,:).*0.5;
            
    end
    scatterSet{i}.MarkerEdgeColor=ax.ColorOrder(mod(i-1,size(ax.ColorOrder,1))+1,:);
    if ~isempty(EdgeColor)
        scatterSet{i}.LineWidth=1.5;
        scatterSet{i}.MarkerEdgeColor=EdgeColor;
        scatterSet{i}.SizeData=60;
        scatterSet{i}.MarkerFaceColor=ax.ColorOrder(mod(i-1,size(ax.ColorOrder,1))+1,:);
    end
end
for i=1:length(lineSet)
    switch theme.EdgeStyle
        case 'none'
            EdgeColor=[];
        case 'gray'
            EdgeColor=[0.3 0.3 0.3];
        case 'white'
            EdgeColor=[0.96 0.96 0.96];
        case 'ori'
            EdgeColor=ax.ColorOrder(mod(i-1,size(ax.ColorOrder,1))+1,:).*0.5;
            
    end
    lineSet{i}.LineWidth=1.8;
    lineSet{i}.Color=ax.ColorOrder(mod(i-1,size(ax.ColorOrder,1))+1,:);
    if ~isempty(EdgeColor)
        lineSet{i}.LineWidth=1.5;
        lineSet{i}.MarkerEdgeColor=EdgeColor;
        lineSet{i}.MarkerSize=8;
        lineSet{i}.MarkerFaceColor=ax.ColorOrder(mod(i-1,size(ax.ColorOrder,1))+1,:);
    end
end
end
    


% legend  Stylized 
if ~isempty(ax.Legend)
switch theme.LegendStyle
    case 'ggplot'
        ax.Legend.FontSize=11;
        ax.Legend.Title.FontSize=14;
        ax.Legend.AutoUpdate='off';
        if ~isempty(regexpi(ax.Legend.Location,'out'))
            ax.Legend.Box='off';
            lgdPos=ax.Legend.Position;
            xyMin=[(lgdPos(1)-ax.Position(1))/ax.Position(3)*(ax.XLim(2)-ax.XLim(1))+ax.XLim(1),...
                (lgdPos(2)-ax.Position(2))/ax.Position(4)*(ax.YLim(2)-ax.YLim(1))+ax.YLim(1)];
            xyMax=[(lgdPos(1)+lgdPos(3)-ax.Position(1))/ax.Position(3)*(ax.XLim(2)-ax.XLim(1))+ax.XLim(1),...
                (lgdPos(2)+lgdPos(4)-ax.Position(2))/ax.Position(4)*(ax.YLim(2)-ax.YLim(1))+ax.YLim(1)];
            ax.Legend.Title.Visible='off';
            xyMin(1),xyMax(2)
            ax.Legend.UserData.NewTitle=text(ax,xyMin(1),xyMax(2),[' ',ax.Legend.Title.String],...
                'FontSize',14,'VerticalAlignment','top','FontWeight','bold');
        else
            ax.Legend.Box='off';
            lgdPos=ax.Legend.Position;
            xyMin=[(lgdPos(1)-ax.Position(1))/ax.Position(3)*(ax.XLim(2)-ax.XLim(1))+ax.XLim(1),...
                (lgdPos(2)-ax.Position(2))/ax.Position(4)*(ax.YLim(2)-ax.YLim(1))+ax.YLim(1)];
            xyMax=[(lgdPos(1)+lgdPos(3)-ax.Position(1))/ax.Position(3)*(ax.XLim(2)-ax.XLim(1))+ax.XLim(1),...
                (lgdPos(2)+lgdPos(4)-ax.Position(2))/ax.Position(4)*(ax.YLim(2)-ax.YLim(1))+ax.YLim(1)];
            xDiff=(xyMax(1)-xyMin(1));
            yDiff=(xyMax(2)-xyMin(2));
            ax.Legend.UserData.NewBkg=rectangle(ax,'Position',[xyMin,xDiff,yDiff],'Curvature',0.2,...
                'LineWidth',1.2,'EdgeColor',[0.39 0.41 0.39],'FaceColor',[1 1 1 .2]);
            %ax.Legend.Title.FontSize=14;
            ax.Legend.Title.Visible='off';
            ax.Legend.UserData.NewTitle=text(ax,xyMin(1),xyMax(2),[' ',ax.Legend.Title.String],...
                'FontSize',14,'VerticalAlignment','top','FontWeight','bold');
        end
    case 'own1'
        ax.Legend.Color=[0.9412 0.9412 0.9412];
        ax.Legend.LineWidth=0.8;
        ax.Legend.FontSize=11;
end
end

% axes Stylized 
switch theme.AxesTheme
    case 'gray'
       ax.Parent.Color=[1 1 1];
       ax.Color=[0.9,0.9,0.9];
       ax.Box='off';
       grid(ax,'on');
       ax.TickDir='out';
       ax.GridColor=[1 1 1];
       ax.GridAlpha=1;
       ax.LineWidth=1.2;
       ax.XColor=[0.33,0.33,0.33];
       ax.YColor=[0.33,0.33,0.33];
       ax.TickLength=[0.015 0.025];
       plot(ax,[ax.XLim(2),ax.XLim(1),ax.XLim(1),ax.XLim(2),ax.XLim(2)],...
               [ax.YLim(2),ax.YLim(2),ax.YLim(1),ax.YLim(1),ax.YLim(2)],...
               'Color',[1 1 1],'LineWidth',2)
    case 'economist'
        ax.Parent.Color=[0.8400 0.8900 0.9200];
        ax.Color=[0.8400 0.8900 0.9200];
        ax.Parent.InvertHardcopy='off';
        ax.Box='off';
        ax.YGrid='on';
        ax.GridColor=[1 1 1];
        ax.GridAlpha=1;
        ax.LineWidth=1.2;
        ax.XColor=[0.33,0.33,0.33];
        ax.YColor='none';
        ax.TickLength=[0.015 0.025];
        for i=1:length(ax.YTick)
            ax.UserData.NewYTick(i)=...
                text(ax,ax.XLim(1)-ax.TickLength(1)/(ax.Position(3))*(ax.XLim(2)-ax.XLim(1)),...
                ax.YTick(i),ax.YTickLabel{i},'HorizontalAlignment','right','Color',[0.33,0.33,0.33]);
        end
    case 'wsj'
        ax.Parent.Color=[0.9700 0.9500 0.8900];
        ax.Color=[0.9700 0.9500 0.8900];
        ax.Parent.InvertHardcopy='off';
        ax.Box='off';
        ax.YGrid='on';
        ax.GridAlpha=1;
        ax.LineWidth=0.8;
        ax.YColor='none';
        ax.GridLineStyle=':';
        ax.TickLength=[0.015 0.025];
        for i=1:length(ax.YTick)
            ax.UserData.NewYTick(i)=...
                text(ax,ax.XLim(1)-ax.TickLength(1)/(ax.Position(3))*(ax.XLim(2)-ax.XLim(1)),...
                ax.YTick(i),ax.YTickLabel{i},'HorizontalAlignment','right','Color',[0.33,0.33,0.33]);
        end
    case 'own1'
        grid(ax,'on');
        ax.GridLineStyle='--';
        ax.LineWidth = 1;
end
end