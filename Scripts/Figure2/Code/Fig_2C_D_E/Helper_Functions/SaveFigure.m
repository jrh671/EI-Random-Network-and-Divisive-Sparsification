figList = findall(groot, 'Type', 'Figure'); % Find all open figures

% Loop through each figure
for i = 1:length(figList)
    fig = figList(i);
    
    % Set your desired output filename
    filename = sprintf('Figure_%d.png', fig.Number);
    
    % Save the figure
    saveas(fig, filename);
end
