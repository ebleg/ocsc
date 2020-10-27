function [ax] = plot_traffic_simulation(fig, u, x0, label, par)
    [~, x] = J(u, x0, par);
    
    t = 1:60;

    axes = flip(findall(fig,'type','axes'));
    
    titles = {'Link $ud$ - Queue for $o_1$', 'Link $o_1d$ - Queue for $o_2$', ...
              'Link $ud$ - Queue for $o_2$' 'Link $o_1d$ - Queue for $o_3$', ...
              '$n_{ud}$' '$n_{o_1d}$'};

    if isempty(axes) % No plots yet
        tiles = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
%       title('\textbf{State evolution over time}')
        for i = 1:6
           axes(i) = nexttile(tiles);
           hold(axes(i), 'on');
        end
        title(tiles, '\textbf{State evolution over time}', 'interpreter', 'latex');
        xlabel(tiles, 'Time [min]', 'interpreter', 'latex');
        ylabel(tiles, 'Number of cars [-]', 'interpreter', 'latex');
    end

    for state = 1:4
        for link = 1:2
            if state ~= 3
                lin_idx = min([(state-1), 2])*2 + link;
                plot(axes(lin_idx), t, x(state + 4*(link-1), :), ...
                    'DisplayName', label, 'LineWidth', 1)
                title(axes(lin_idx), titles{lin_idx});
            end
        end
    end
    
    % Need Matlab R2020b
    lgd = legend(axes(1));
    lgd.Orientation = 'horizontal';
    lgd.Layout.Tile = 'north';
    
end

