function downstream_neuron_output = RRN_DownstreamNeuron(spikeMat,TopNeurons)
 
ChosenNeurons = TopNeurons;

% Group of E neurons of interest
group = ChosenNeurons;

% Define a percentage threshold for co-activity (e.g., 60%)
percentage_threshold = .6;

% Define firing thresholds
FRLow = 0;
FRHigh = 1;

% The number of neurons in the group that must co-fire for the downstream neuron to fire
coactivity_threshold = round(length(group) * percentage_threshold);

% Create a variable to store the output of the downstream neuron
downstream_neuron_output = 0; % initialize to 0

% Calculate the total firing rate for each neuron in the group over the 5-second interval
firing_rates = sum(spikeMat(group, :), 2) / size(spikeMat, 2);

% Check if the number of neurons firing above the rate threshold exceeds the coactivity threshold
if sum(firing_rates > FRLow) >= coactivity_threshold && sum(firing_rates < FRHigh)
    downstream_neuron_output = 1;
end
end