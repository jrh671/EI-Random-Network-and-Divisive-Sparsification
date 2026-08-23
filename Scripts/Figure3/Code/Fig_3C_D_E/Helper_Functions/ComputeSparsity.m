
%% Fast Tm Sparsity

[Rate,order]=Compute_RateMap(FiringRateF,Positions);

zero_mask = Rate == 0;

sparsity_per_time = sum(zero_mask, 1) / size(Rate, 1);

SparsityFiringF = mean(sparsity_per_time);

%% Slow Tm Sparsity

[Rate,order]=Compute_RateMap(FiringRateS,Positions);

zero_mask = Rate == 0;

sparsity_per_time = sum(zero_mask, 1) / size(Rate, 1);

SparsityFiringS = mean(sparsity_per_time);

%% DivSparse Sparsity
[Active,order]=Compute_RateMap(Activity,Positions);

% Step 1: Create a binary mask where elements are 1 if the value is 0
zero_mask = Active == 0;

% Step 2: Compute the proportion of rows that are 0 for each column
sparsity_per_time = sum(zero_mask, 1) / size(Active, 1);

% Step 3: Compute the average sparsity across all columns
SparsityActivity = mean(sparsity_per_time);