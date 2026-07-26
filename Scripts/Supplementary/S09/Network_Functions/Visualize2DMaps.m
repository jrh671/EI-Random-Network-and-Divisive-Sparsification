
    firing_rate_snapshots_CA3 = FiringRates;
    if VisualizePlasticity==0
    ComputeMaps; 
    end

    if VisualizePositions == 1
    
    figure;
    hold on;
    
    % Normalize time indices to [0,1] for colormap scaling
    num_points = size(trajectory, 1);
    colors = spring(num_points); % Get colors from the 'spring' colormap
    
    % Plot each segment individually with its corresponding color
    for i = 1:num_points-1
        plot(trajectory(i:i+1, 1), trajectory(i:i+1, 2), '-', 'Color', colors(i, :), 'LineWidth', 1.5);
    end
    
    SheetPOS=[trajectory(:,1),trajectory(:,2)];

    hold off;
    colormap(spring); % Apply colormap
    colorbar; % Show color scale
    caxis([0 1]); % Ensure colormap range is 0 to 1
    
    title('Trajectory Colored by Time');
    xlabel('X Position');
    ylabel('Y Position');

        pause(1)
    end
  
    if VisualizePlasticity==1
        chosen_neurons=Choose2NeuronsCA3;
        % Visualization of the chosen neurons' activities evolving through time
        figure('units','normalized','outerposition',[0 0 1 1])

         for t = StartTime:EndTime
             Neurons=chosen_neurons;
             for n = 1:2
                    subplot(1, 2, n);
                    p_t = Traj1(t, :);

                    % chosen_neuron = chosen_neurons_CA3(end-2+n);
                    chosen_neuron=Neurons(n);
                    activity = firing_rate_snapshots_CA3{t}{chosen_neuron};
                    imagesc(activity);

                    set(gca, 'YDir', 'normal');  % Make y-axis upright
                    axis equal
                    axis tight
                    hold on;
                    plot(p_t(2), p_t(1), 'rx', 'MarkerSize', 25, 'LineWidth', 2);  % Mark the current position
                    if special_radius > 0
                        plot(special_location(2), special_location(1), 'go', 'MarkerSize', 8, 'LineWidth', 2);  % Mark the special location
                        viscircles(fliplr(special_location), special_radius, 'Color', 'g'); % Draw the special location radius
                    end
                    title(sprintf('Neuron %d Activity at Time %d', chosen_neuron, t));
                    xlim([1, P])
                    ylim([1, P])
                    clim([0,10])

                    hold off;
                    colorbar;

            end
            pause(0.1);  % Pause to create an animation effect
         end

              times=[1,30,60,90,100];

          for t = 1:5
              figure
             Neurons=Choose2NeuronsCA3;
             Time=times(t);
             for n = 1:2
                    subplot(1, 2, n);
                        p_t = Traj1(Time, :);
    
                    % chosen_neuron = chosen_neurons_CA3(end-2+n);
                    chosen_neuron=Neurons(n);
                    activity = firing_rate_snapshots_CA3{Time}{chosen_neuron};
                    SheetRates{t,n}=activity;
                    imagesc(activity);
                    set(gca, 'YDir', 'normal');  % Make y-axis upright
                    axis equal
                    axis tight
                    hold on;
                    plot(p_t(2), p_t(1), 'rx', 'MarkerSize', 25, 'LineWidth', 2);  % Mark the current position
                    if special_radius > 0
                        plot(special_location(2), special_location(1), 'go', 'MarkerSize', 8, 'LineWidth', 2);  % Mark the special location
                        viscircles(fliplr(special_location), special_radius, 'Color', 'g'); % Draw the special location radius
                    end
                    title(sprintf('Neuron %d Activity at Time %d', chosen_neuron, Time));
                    xlim([1, P])
                    ylim([1, P])
                    clim([0,10])
                    hold off;
                    colorbar;

            end
            pause(0.5);  % Pause to create an animation effect
            clim([0,10])

          end
          
    end

MCells=find(memory_neurons_CA3==1);

