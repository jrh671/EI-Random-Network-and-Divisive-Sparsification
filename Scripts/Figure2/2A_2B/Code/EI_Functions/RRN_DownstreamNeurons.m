function downstream_neuron_outputs = RRN_DownstreamNeurons(spikeMat,TopNeurons)

neuron_indices = TopNeurons;

% Define a percentage threshold for co-activity (e.g., 90% --> 0.90)
percentage_threshold = .7;

% Define a firing rate threshold (e.g., 0.1, meaning neuron should fire at least 10% of the time)
FR = 0;

% The number of neurons in the group that must co-fire for the downstream neuron to fire
coactivity_threshold = round(20 * percentage_threshold);

% Create a variable to store the output of the downstream neuron
downstream_neuron_outputs = zeros(25,1); % initialize to 0

% Loop over the 25 downstream neurons
for i = 1:25
    % Define the group of E neurons of interest for the current downstream neuron
    group = neuron_indices((i-1)*20+1 : i*20);

    % Calculate the total firing rate for each neuron in the group over the 5-second interval
    firing_rates = sum(spikeMat(group, :), 2) / size(spikeMat, 2);

    % Check if the number of neurons firing above the rate threshold exceeds the coactivity threshold
    if sum(firing_rates > FR) >= coactivity_threshold
        downstream_neuron_outputs(i) = 1;
    end
end




end