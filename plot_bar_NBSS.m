
function plot_bar_NBSS(ax,data,colors, app)
cla(ax)
% figure
% ax=axes;


% data=base_Zooscan(1).regroupped.Ybv_Ellipsoid_BV_spectra;
% X=exp(base_Zooscan(1).tot.X)  ;
% ESD=2*(X*3/(4*pi)).^(1/3); % in mm
%     ESD=ESD*1000; %in Âµm
%     
%     size(log10(app.ESDsize))
%     size(data)

switch app.NBSSSwitch.Value
    case 'surface'
    
    h=area(ax,app.ESDsize,data,'FaceColor','flat');
    for k = 1:size(data,2)
    h(k).CData = k;
end
    ax.XScale='log';
%     colors = (jet(size(data,2))+1)/2;
% colors =linspecer(169,'qualitative')
% colors=[vega20(16)]
ax.XLim=[app.ESDminEditField.Value app.ESDmaxEditField.Value];
colormap(ax,colors)
ax.XLabel.String='ESD (µm)';
ax.YLabel.String='NBSS mm^{3} mm^{-3} m^{-3}';

    case 'bar'

h=bar(ax,log10(app.ESDsize),data,'stacked','FaceColor','flat')  ;
minors=[0.2:0.1:0.9 2:9 20:10:90 200:100:900 2000:1000:9000 20000:10000:90000];

ax.XTick=[-1:6]; %// adjust manually; values in log scale
ax.XTickLabel=10.^(ax.XTick); %// use labels with linear values
%set(gca,'XMinorTick','on','XMinorTickValues',minors)
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues=log10(minors);
for k = 1:size(data,2)
    h(k).CData = k;
end
% colormap(linspecer(169))
%set(h,'Edgecolor',[0.9 0.9 0.9])
ax.XLim=[log10(app.ESDminEditField.Value) log10(app.ESDmaxEditField.Value)];
% colors = (jet(size(data,2))+1)/2;
% colors =linspecer(169,'qualitative')
% colors=[vega20(16)]
colormap(ax,colors)
ax.XLabel.String='ESD (µm)';
ax.YLabel.String='NBSS mm^{3} mm^{-3} m^{-3}';

end
