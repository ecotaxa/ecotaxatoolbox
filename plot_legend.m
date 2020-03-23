function plot_legend(ax,data,list,colors,app)

[n,m]=size(data);
data=ones(n,m);

% figure
% ax=axes
cla(ax)
colorbar(ax,'off');
%             hold(ax,'on')
%             h=pie(ax,100*data/sum(data));
%             
% %             
% %             h=bar(ax,data','stacked','FaceColor','flat')  ;

set(ax,'Clim',[1 m])
%n=length(data);
%colors = (jet(n)+1)/2;
                            %colors=[vega20b(20);vega20c(20)];
                            colormap(ax,colors) 
                            h=colorbar(ax);
                            set(ax,'visible','off')
                            
                           get(h,'position')
                           if m>25
                               fontcolormap=8;
                           else
                               fontcolormap=12;
                           end
                     set(h,'ticks',1:m,'ticklabels',list,'fontsize',fontcolormap,'position',[0.1 0.05 0.1 0.9],'AxisLocation','in')       
%                             legendhandle=h(1:2:end);
%     legendnames=list;
%     ah1 = gca;
%     % Legend at axes 1
%     legend(ax,legendhandle(1:50),legendnames(1:50),'location','westoutside');
% 
%     % Block 2
%     % Axes handle 2 (unvisible, only for place the second legend)
%     ah2=axes('position',get(gca,'position'), 'visible','off');
%     % Legend at axes 2
%     legend(ah2,legendhandle(51:n),legendnames(51:n),'location','eastoutside','NumColumns',2);
% 
%     
    
    
%legend(ax,list)
% for k = 1:size(h,2)
%     h(k).Visible = 'off';
% end


