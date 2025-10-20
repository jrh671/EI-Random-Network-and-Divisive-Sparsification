if VisualizePlasticity==1
PlottingEpochs=2;
else
PlottingEpochs=Epochs;
end

for Epoch = PlottingEpochs
    sorted_indices_CA3 = Sorting{Epoch};
    firing_rate_snapshots_CA3 = Rates{Epoch};
    
    if Epoch == 2 || Epoch ==4
        special_radius = GoalRadius; 
    else
        special_radius = 0;
    end
    % Find the top 4 neurons with the highest average activity at the special location for CA1
    chosen_neurons_CA3 = sorted_indices_CA3(1:8);
    chosen_neurons_CA3(length(chosen_neurons_CA3)-1:length(chosen_neurons_CA3)) = Choose2NeuronsCA3;
    
    ComputeMaps; 

    ActiveRates = temporal_rates_CA3 > 3*mu_rho/4;
    TemporalRates{Epoch,1}=temporal_rates_CA3;
    
    if VisualizePositions == 1
    
        figure;
        plot(trajectory(:, 1), trajectory(:, 2), 'o-');
        if special_radius > 0
            hold on;
            plot(special_location(1), special_location(2), 'rx', 'MarkerSize', 12, 'LineWidth', 2);  % Mark the special location
            viscircles(special_location, special_radius, 'Color', 'r'); % Draw the special location radius
            hold off;
        end
        xlabel('X Position');
        ylabel('Y Position');
        title('Full Combined Trajectory');
        grid on;
        axis equal
        axis tight
        pause(1)
    end
  
    if VisualizePlasticity==1
        chosen_neurons=Choose2NeuronsCA3;
        if Epoch == 2 || Epoch == 4 || Epoch == 6 || Epoch == 8
        % Visualization of the chosen neurons' activities evolving through time
        figure('units','normalized','outerposition',[0 0 1 1])

         for t = StartTime:EndTime
             Neurons=chosen_neurons;
             for n = 1:2
                    subplot(1, 2, n);
                    if Context==1
                        p_t = Trajectory1(t, :);
                    elseif Context==2
                        p_t = Trajectory2(t, :);
                    end
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


          for t = [1,11,29,50,100]
              figure
             Neurons=Choose2NeuronsCA3;
             for n = 1:2
                    subplot(1, 2, n);
                    if Context==1
                        p_t = Trajectory1(t, :);
                    elseif Context==2
                        p_t = Trajectory2(t, :);
                    end
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
            pause(0.5);  % Pause to create an animation effect
            clim([0,10])

          end
          
        end
    end
end

MCells=find(memory_neurons_CA3==1);


if VisualizePlasticity==0
if Context==1
    Epoch = 4;
Anti_Cofiring;
SheetCtx1Epoch1=sortedKendallTauMatrix;
Sort=sortedIndices;
Memory=activeMemoryCells;
   Epoch = 6;

Anti_Cofiring;

    SheetCtx1Epoch2=sortedKendallTauMatrix;

    Epoch = 2;
elseif Context==2
    Epoch = 4;

Anti_Cofiring;
SheetCtx2Epoch2=sortedKendallTauMatrix;

Sort=sortedIndices;
Memory=activeMemoryCells;
    Epoch = 6;
Anti_Cofiring;
SheetCtx2Epoch1=sortedKendallTauMatrix;

    Epoch = 2;
end

MapStability;
end

