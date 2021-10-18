function plot_fit_truncated(x,y,fitted_y,x_truncated, fitted_y_truncated,lim, trace)
        figure1=figure( 'Name', ' fits' );
        % Create axes
        axes1 = axes('Parent',figure1);
        hold(axes1,'on');
        
        plot(x(1:length(x),1),y(1:length(y),1),'DisplayName','avgFRET vs. time(sec)','LineWidth',0.5,...
        'Color',[0 0 1]);
        plot(x(1:length(x),1),fitted_y(1:length(fitted_y),1),'DisplayName','ktot(full turnover) fit','LineWidth',1.5,...
        'Color',[1 0 0]);
    
        plot(x_truncated,fitted_y_truncated,'g*');
        title({'trace',trace},'FontSize',24);
       % Create xlabel
        xlabel('Time(sec)');
        % Create ylabel
        ylabel('FRET','Interpreter','none');
        xlim(axes1,[0 lim]);
        ylim(axes1,[0.4 1]);
        box(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',24);
        % put waiting time and radius of liposome to plot
       %ylim_text=get(gca,'ylim');
       %xlim_text=get(gca,'xlim');
       %txt=sprintf('liposome size = %.1f nm\ntransport time = %.1f sec', r,vmax);
       %text(xlim_text(2)-490,ylim_text(2)-0.08,txt,'FontSize',14);
       fname='fits';
       filename=sprintf('trace %d',trace);
       saveas(figure1,fullfile(fname, filename),'png')

end