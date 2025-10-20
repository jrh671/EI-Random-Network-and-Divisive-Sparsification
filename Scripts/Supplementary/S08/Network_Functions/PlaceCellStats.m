% Place Cells
meanValues = cellfun(@(x) mean(x(:)), firing_rates_CA3);
meanValues(meanValues==10)=0;
PlaceCells=find(meanValues>0);
numPlaceCells = N-sum(meanValues==0);

% Assuming your cell array is called 'cellArray'
numPeaks = zeros(N, 1); % Preallocate for storing peak counts

sigma = 2; % Set smoothing intensity
for i = 1:N
    CA{i} = imgaussfilt(firing_rates_CA3{i}, sigma);
end

for i = 1:N
    matrix = CA{i}; % Extract the 7x7 matrix
    
    % Compute the gradient magnitude
    [Gx, Gy] = gradient(matrix); 
    gradMagnitude = sqrt(Gx.^2 + Gy.^2);
    
    % Identify peaks using regional maxima
    peaks = imregionalmax(matrix); 
    
    % Count the number of peaks
    numPeaks(i) = sum(peaks(:));
end

numPeaks(numPeaks>10)=0;

TwoPeaks = find(numPeaks==2);
length(TwoPeaks)
figure;
histogram(numPeaks, 'BinEdges', 0:1:4, 'EdgeColor', 'black', 'FaceColor', 'blue');
% xlim([0, 4]);
ylim([0, 80]);
% xticks(0.5:1:3.5); % Set tick positions at bin centers
% xticklabels({'0', '1', '2', '3'}); % Label them as integer values
title('Place Field Distribution');
xlabel('Number of PFs');
ylabel('Neuron Count');
