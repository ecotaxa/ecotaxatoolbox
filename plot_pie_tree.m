function plot_pie_tree(ax,data,list,colors,app)


cla(ax)
            hold(ax,'on')
            
            switch app.proportionsSwitch.Value
                
                case 'Pie'
                    axis(ax,'normal')
                    if app.TextCheckBox.Value==0
                        xlim(ax,'auto');
                        falsetext={' '};
                               n=length(data);
                               falsetext=repmat(falsetext,[n,1]);
                            pie(ax,100*data/sum(data),falsetext);
                            set(ax,'xlim',[-1 1],'ylim',[-1 1])

                    else
                           if app.PercentagesCheckBox.Value==1

                               pie(ax,100*data/sum(data));%app.listspeciesshort); %,app.listspeciesshort,
                           set(ax,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5])
                           else
%                                falsetext={' '};
%                                n=length(data);
%                                falsetext=repmat(falsetext,[n,1]);
                            h=pie(ax,100*data/sum(data),list); %,app.listspeciesshort,
                            %delete(h(2:2:end)) %deleting labels
                            set(ax,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5])
                           end
                    end
%                     n=length(data);
%                             colors = (jet(n)+1)/2;
                            %colors=[vega20b(20);vega20c(20)];
                            colormap(ax,colors)
                            set(ax,'visible','off')
                                axis(ax,'equal')
    %axis(ax,'tight')
    %axis(ax,'off')
                    
                    %legend(ax,list)
                    
                case 'Treemap'
                    
%                     n=length(data);
%                     colors = (jet(n)+1)/2;
                    rectangles = treemapui(data);
                    if app.TextCheckBox.Value==0
                        
                            plotRectanglesui(ax,rectangles,[],colors);
                            
                    else
                            plotRectanglesui(ax,rectangles,list,colors);
                    end
                    
            end
            

            