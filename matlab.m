% Student License -- for use by students to meet course requirements
% and perform academic research at degree granting institutions only.

% Use your loaded table
T = readtable('C:\Users\walsh\Downloads\CRYO-SEM Accuracy INTERNAL\CRYO-SEM X10000\CRYO-SEM X10000 [1]\CRYO-SEM X10000 [1] STATS\GOLD STANDARD [X10000] Results.csv');

% Part 1: Calculate Euclidean Distances

% Use capital 'X' and 'Y' from your table
X = T.X;
Y = T.Y;
n = length(X);

% Initialize distance matrix
DistanceMatrix = NaN(n, n);

% Inline Euclidean distance formula — no function needed
for i = 1:n
    for j = 1:n
        dx = X(i) - X(j);
        dy = Y(i) - Y(j);
        DistanceMatrix(i, j) = sqrt(dx^2 + dy^2);
    end
end

% Create labeled table
rowNames = strcat("Pore_", string(1:n));
DistanceTable = array2table(DistanceMatrix, ...
    'VariableNames', rowNames, 'RowNames', rowNames);

% Save distance matrix to GitHub folder
distOutPath = 'C:\Users\walsh\Documents\GitHub\AGAROSE-HYDROGEL-TRENDS-USING-AI-ML\M_1000xc_euclidean_distances.csv';
writetable(DistanceTable, distOutPath, 'WriteRowNames', true);
fprintf('Saved: %s\n', distOutPath);

% Extract upper triangle (excluding diagonal)
dists = [];
for i = 1:n
    for j = i+1:n
        dists(end+1) = DistanceMatrix(i, j); %#ok<SAGROW>
    end
end

% Summary Statistics 

% Count pores
poreCount = height(T);
fprintf('Total pore count: %d\n', poreCount);

% AECD (Area-Equivalent Circular Diameter)
AECD = sqrt(4 * T.Area / pi);

% Exclude variables from summary
excludeVars = {'x', 'X', 'Y', 'FeretX', 'FeretY', 'FeretAngle', 'Angle', 'AECD'};
allVars = T.Properties.VariableNames;

metrics = {};
meanVals = [];
ciLower = [];
ciUpper = [];

% 95% confidence interval function
computeCI = @(data) deal( ...
    mean(data, 'omitnan') - tinv(0.975, sum(~isnan(data)) - 1) * std(data, 'omitnan') / sqrt(sum(~isnan(data))), ...
    mean(data, 'omitnan') + tinv(0.975, sum(~isnan(data)) - 1) * std(data, 'omitnan') / sqrt(sum(~isnan(data))) );

% Loop through all variables
for i = 1:numel(allVars)
    varName = allVars{i};
    if isnumeric(T.(varName)) && ~ismember(varName, excludeVars)
        val = T.(varName);
        metrics{end+1} = varName; %#ok<SAGROW>
        meanVals(end+1) = mean(val, 'omitnan'); %#ok<SAGROW>
        [low, high] = computeCI(val);
        ciLower(end+1) = low; %#ok<SAGROW>
        ciUpper(end+1) = high; %#ok<SAGROW>
    end
end

% Add AECD
metrics{end+1} = 'AECD';
meanVals(end+1) = mean(AECD, 'omitnan');
[low, high] = computeCI(AECD);
ciLower(end+1) = low;
ciUpper(end+1) = high;

% Add Euclidean distance
metrics{end+1} = 'EuclideanDistance';
meanVals(end+1) = mean(dists, 'omitnan');
[low, high] = computeCI(dists);
ciLower(end+1) = low;
ciUpper(end+1) = high;

% Round values
meanVals = round(meanVals, 4);
ciLower = round(ciLower, 4);
ciUpper = round(ciUpper, 4);

% Save summary
summaryTable = table(metrics', meanVals', ciLower', ciUpper', ...
    'VariableNames', {'Metric', 'Mean', 'CI95_Lower', 'CI95_Upper'});
disp('--- Final Summary ---');
disp(summaryTable);

% Save summary to GitHub folder
statsOutPath = 'C:\Users\walsh\Documents\GitHub\AGAROSE-HYDROGEL-TRENDS-USING-AI-ML\M_1000xc_stats.csv';
writetable(summaryTable, statsOutPath);
fprintf('Saved: %s\n', statsOutPath);
