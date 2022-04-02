function plot_projections(projections, H, T, start_color, end_color, ylabel)
    ax = gca;
    fig.Units = "centimeters";
    fig.Position = [0 0 18 4];
    fig.PaperUnits = "centimeters";
    fig.PaperPosition = [0 0 18 4];
    fig.PaperSize = [18 4];
    ax.Units = "centimeters";
    ax.Position = [0 0 20 4];
    ax.OuterPosition = [-1 0 20 4];
    hold on;
    color = colorGradient(start_color, end_color, H);
    for h=2:H
        plot(projections(h,:), "Color", color(h,:), 'LineStyle', ':');
    end
    plot(projections(1,:), "Color", color(1,:), "LineWidth", 2);
    for t=1:T
        xline(24*t, "Color", "k", "LineStyle", "--");
        text(24*t - 0.5, max(projections,[],"all") - 0.5, compose('Day %d', t), ...
             "FontSize", 8, "HorizontalAlignment", "right");
    end
    ax.XLimitMethod = "tight";
    ax.YLabel.String = ylabel;
    grid on;
end

