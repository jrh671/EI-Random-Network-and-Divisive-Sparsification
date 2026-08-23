
SaveMean=cell(8,1);

% 1: Assembly Tagging | 2: Template | 3: L1| 4: L2| 5: Linear

for Decoder = 1 %Out of 5
    if ShowSparsity == 0 && FirstPlot==1
figure;
    end
Instance=1;
RunDecoder_Stats_3;
end
