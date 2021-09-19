%% Load GM data and convert to matrix
load('public_gm12878_reads.mat')
binSize = 1e5; % 100kb
public_gm12878_reads.posA = ceil(public_gm12878_reads.posA/binSize);
public_gm12878_reads.posB = ceil(public_gm12878_reads.posB/binSize);
public_gm12878_strands = sortrows(public_gm12878_reads, 'read_id', 'ascend');
porecMatrix = table2array(public_gm12878_reads);

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
porecMatrix = porecMatrix(uniqReadContacts, :);
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

incidenceMatrix = incidenceMatrix(115:129, :); % Select the region between 115 and 129
incidenceMatrix = unique(incidenceMatrix', 'rows')'; % Remove duplicate hyperedges

%% Extract all higher-order hyperedges
highOrderContacts = cell(size(incidenceMatrix, 2), 1);
for i = 1:size(incidenceMatrix, 2)
    highOrderContacts{i} = find(incidenceMatrix(:, i)>0);    
end
highOrderContacts = flipud(highOrderContacts(cellfun('length', highOrderContacts)>=3)); % Order >=3
highOrderContacts = cellfun(@(v) v+114, highOrderContacts, 'UniformOutput', false); % Add the offset 114

%% Create CSV file for PAOHvis
CSVIncidenceMatrix = [];
label = [];
for i = 1:length(highOrderContacts)
    contact = highOrderContacts{i};
    for j = 1:length(contact)
        CSVIncidenceMatrix = [CSVIncidenceMatrix; [i contact(j)]]; %#ok<AGROW>
        if contact(j) >= 118 && contact(j) <= 123
            label{end+1} = 'TAD 1'; %#ok<*SAGROW>
        elseif contact(j) >= 124 && contact(j) <= 127
            label{end+1} = 'TAD 2';
        else 
            label{end+1} = 'Other';
        end
    end 
end 

CSVIncidenceMatrixTable = table(CSVIncidenceMatrix(:, 1), num2str(CSVIncidenceMatrix(:, 2), '%03d'), ...
    cell(size(CSVIncidenceMatrix, 1), 1), cell(size(CSVIncidenceMatrix, 1), 1), label');
writetable(CSVIncidenceMatrixTable, 'Intra_Inter_TAD_GM.csv')
