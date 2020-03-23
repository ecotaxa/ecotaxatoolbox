        function refreshplot(app)
        
        epsilon=10^(-50);
        
       plot_pie_tree(app.UIAxes,app.subsetdataAbshort+epsilon,app.listspeciesshort,app.colormapspecies,app)
              
                
        plot_pie_tree(app.UIAxes_4,app.subsetdataBvshort+epsilon,app.listspeciesshort,app.colormapspecies,app)
            
       plot_pie_tree(app.UIAxes_2,app.subsetdataAbpftshort+epsilon,app.listpftshort,app.colormappft,app)
        
        plot_pie_tree(app.UIAxes_5,app.subsetdataBvpftshort+epsilon,app.listpftshort,app.colormappft,app)
        
       plot_pie_tree(app.UIAxes_3,app.subsetdataAbtrophicshort+epsilon,app.listtrophicshort,app.colormapspecies,app)
        
        plot_pie_tree(app.UIAxes_6,app.subsetdataBvtrophicshort+epsilon,app.listtrophicshort,app.colormapspecies,app)
        
        if app.NBSSCheckBox.Value==1
        
        plot_bar_NBSS(app.UIAxes_7,app.subsetdataNBSSshort,app.colormapspecies,app)
        plot_bar_NBSS(app.UIAxes_8,app.subsetdataNBSSpftshort,app.colormappft,app)
       plot_bar_NBSS(app.UIAxes_9,app.subsetdataNBSStrophicshort,app.colormapspecies,app)
        else
            cla(app.UIAxes_7)
            cla(app.UIAxes_8)
            cla(app.UIAxes_9)
        end
        

          plot_legend(app.UIAxes2,app.subsetdataAbshort+epsilon,app.listspeciesshort,app.colormapspecies,app)
          plot_legend(app.UIAxes2_2,app.subsetdataAbpftshort,app.listpftshort,app.colormappft,app)