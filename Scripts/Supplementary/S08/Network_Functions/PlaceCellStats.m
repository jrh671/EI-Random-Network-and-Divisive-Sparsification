% Place Cells
maxValues = cellfun(@(x) max(x(:)), M_prim_CA3);
maxValues(maxValues==10)=0;
PlaceCells=find(maxValues>0);
numPlaceCells = N-sum(maxValues==0);

% Assuming your cell array is called 'cellArray'
numPeaks = zeros(N, 1); % Preallocate for storing peak counts
% 
% sigma = 0; % Set smoothing intensity
% for i = 1:N
%     CA{i} = imgaussfilt(M_prim_CA3{i}, sigma);
% end

CA=M_prim_CA3;
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
ThreePeaks = find(numPeaks==3);

figure;
numPeaks(numPeaks == 0) = NaN;
histogram(numPeaks, 'BinEdges', 0:1:14, 'EdgeColor', 'black', 'FaceColor', 'blue');

if Small==1
SheetPeaksSmall=numPeaks;
else
SheetPeaksLarge=numPeaks;
end

% xlim([0, 4]);
ylim([0, 50]);
% xticks(0.5:1:3.5); % Set tick positions at bin centers
% xticklabels({'0', '1', '2', '3'}); % Label them as integer values
title('Place Field Distribution');
xlabel('Number of PFs');
ylabel('Neuron Count');
find(numPeaks>0);
