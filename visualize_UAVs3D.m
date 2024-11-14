%% Visualizing the results
function visualize_UAVs3D(UAV_positions, noof_moves, noof_UAVs, rcom)
    figure;
    % initialize variables needed to draw a sphere 
    [unit_sph_x, unit_sph_y, unit_sph_z] = sphere; 
    
    min_x = min(min(UAV_positions(:,1,:))) - rcom;
    max_x = max(max(UAV_positions(:,1,:))) + rcom;
    min_y = min(min(UAV_positions(:,2,:))) - rcom;
    max_y = max(max(UAV_positions(:,2,:))) + rcom;
    min_z = 0;
    max_z = max(max(UAV_positions(:,3,:))) + rcom;
    x_span = (max_x - min_x);
    y_span = (max_y - min_y);
    x_axis_min = (round(min_x/rcom)-1)*rcom;
    x_axis_max = (round(max_x/rcom)+1)*rcom;
    y_axis_min = (round(min_y/rcom)-1)*rcom;
    y_axis_max = (round(max_y/rcom)+1)*rcom;
    z_axis_min = min_z;
    z_axis_max = (round(max_z/rcom)+1)*rcom;
    iconsize = [x_span/30 y_span/30]; 
    % read image file for UAV icon
    %img_name = 'UAV_icon.png';
    img_name = 'UAV_icon2D.png';
    UAV_icon = imread(img_name, 'png');
    
    for t = 1:noof_moves
       pause(0.25); % pause shortly for viewing
       hold off
       
       for UAV_num = 1:noof_UAVs
            % get the x and y position of UAV
            x_pos = UAV_positions(UAV_num,1,t);
            y_pos = UAV_positions(UAV_num,2,t);
            z_pos = UAV_positions(UAV_num,3,t);

            sph_x = (rcom*unit_sph_x + x_pos);
            sph_y = (rcom*unit_sph_y + y_pos);
            sph_z = (rcom*unit_sph_z + z_pos);
            surf(sph_x, sph_y, sph_z,'FaceColor', 'blue', ...
                 'LineStyle', '-', 'EdgeColor', 'blue', ... 
                 'FaceAlpha', 0.05, 'EdgeAlpha', 0.1);
            xlim([x_axis_min x_axis_max])
            ylim([y_axis_min y_axis_max])
            zlim([z_axis_min z_axis_max])
            hold on
            % put the UAV icon at the UAV's position and adjust opacity
            scatter3(x_pos, y_pos, z_pos, 'filled');
       end
       plot_title = sprintf('UAV positions at t = %d', t);
       title(plot_title);
    end
end