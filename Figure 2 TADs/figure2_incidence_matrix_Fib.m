%% Load Fib data and convert to matrix
load('v1234_np_uniqueIDs.mat')
v1234_np_uniqueIDs = v1234_np_uniqueIDs(:, 2:end);
binSize = 1e5; % 100kb
v1234_np_uniqueIDs.posA = ceil(v1234_np_uniqueIDs.posA/binSize);
v1234_np_uniqueIDs.posB = ceil(v1234_np_uniqueIDs.posB/binSize);
v1234_np_uniqueIDs = sortrows(v1234_np_uniqueIDs, 'read_id', 'ascend');
porecMatrix = table2array(v1234_np_uniqueIDs);

%% Add chromosome offsets
chrLength = [2490 2422 1983 1903 1816 1709 1594 1452 1384 1338 1351 1333 ...
    1144 1071 1020 904 833 804 587 645 468 509 1561 573];
for chr = 2:24
    binAdd = sum(chrLength(1:chr-1));
    binIdxL = porecMatrix(:, 2)==chr;
    porecMatrix(binIdxL, 3) = porecMatrix(binIdxL, 3)+binAdd; 
    binIdxR = porecMatrix(:, 4)==chr;
    porecMatrix(binIdxR, 5) = porecMatrix(binIdxR, 5)+binAdd;
end

%% Threshold pairwise contacts removing noise
[~, uniqReadContacts, ~] = unique(porecMatrix(:, [1 3 5]), 'rows');
porecMatrix_2 = porecMatrix(uniqReadContacts, :);
pairedContacts = porecMatrix(:, [3 5]);
[~, ~, contactIdx] = unique(pairedContacts, 'rows', 'stable');  
contactHist = accumarray(contactIdx, 1);
contactCount = contactHist(contactIdx);
eps = prctile(contactCount, 85); % 85th percentile 
porecMatrix = porecMatrix(contactCount>=eps, :);

%% Extract a chromosome of interest
chr = 22; % Choose chromosome 22
chrIdx = logical((porecMatrix(:, 2)==chr).*(porecMatrix(:, 4)==chr));
porecMatrixChr = porecMatrix(chrIdx, :);
porecMatrixChr(:, [3 5]) = porecMatrixChr(:, [3 5])-sum(chrLength(1:chr-1));

%% Create an incidence matrix
[~, ~, readIDReIdx] = unique(porecMatrixChr(:, 1)); % Reindex the ReadID
porecMatrixChr(:, 1) = readIDReIdx;
incidenceMatrixA = sparse(porecMatrixChr(:, 3), porecMatrixChr(:, 1), 1,...
    chrLength(chr), porecMatrixChr(end, 1));
incidenceMatrixB = sparse(porecMatrixChr(:, 5), porecMatrixChr(:, 1), 1,...
    chrLength(chr), porecMatrixChr(end, 1));
incidenceMatrix = (incidenceMatrixA+incidenceMatrixB)>0;

incidenceMatrix = incidenceMatrix(331:348, :); % Select the region between 115 and 129
incidenceMatrix = unique(incidenceMatrix', 'rows')'; % Remove duplicate hyperedges

%% Extract all higher-order hyperedges
highOrderContacts = cell(size(incidenceMatrix, 2), 1);
for i = 1:size(incidenceMatrix, 2)
    highOrderContacts{i} = find(incidenceMatrix(:, i)>0);    
end
highOrderContacts = flipud(highOrderContacts(cellfun('length', highOrderContacts)>=3)); % Order >=3
highOrderContacts = cellfun(@(v) v+330, highOrderContacts, 'UniformOutput', false); % Add the offset 114

%% Create csv file for PAOHvis
CSVIncidenceMatrix = [];
label = [];
for i = 1:length(highOrderContacts)
    contact = highOrderContacts{i};
    for j = 1:length(contact)
        CSVIncidenceMatrix = [CSVIncidenceMatrix; [i contact(j)]]; %#ok<AGROW>
        if contact(j) >= 334 && contact(j) <= 341
            label{end+1} = 'TAD 1'; %#ok<*SAGROW>
        elseif contact(j) >= 342 && contact(j) <= 346
            label{end+1} = 'TAD 2';
        else 
            label{end+1} = 'Other';
        end
    end 
end 
        
CSVIncidenceMatrixTable = table(CSVIncidenceMatrix(:, 1), num2str(CSVIncidenceMatrix(:, 2), '%03d'), ...
    cell(size(CSVIncidenceMatrix, 1), 1), cell(size(CSVIncidenceMatrix, 1), 1), label');
writetable(CSVIncidenceMatrixTable, 'Intra_Inter_TAD_Fib.csv')
