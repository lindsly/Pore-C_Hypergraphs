%% Load dataset and convert to matrix
load('NlaIII_digest_IRFibroblast_006_reads.mat') % Fibroblast experiment V2 pairs table
binSize = 1e6;
NlaIII_digest_IRFibrobalst_006_reads.posA = ceil(NlaIII_digest_IRFibrobalst_006_reads.posA/binSize);
NlaIII_digest_IRFibrobalst_006_reads.posB = ceil(NlaIII_digest_IRFibrobalst_006_reads.posB/binSize);
NlaIII_digest_IRFibrobalst_006_reads = sortrows(NlaIII_digest_IRFibrobalst_006_reads, 'read_id', 'ascend');
porecMatrix = table2array(NlaIII_digest_IRFibrobalst_006_reads);

%% Add chromosome offsets
chrLength = [249, 243, 199, 191, 182, 171, 160, 146, 139, 134, 136,...
    134, 115, 108, 102, 91, 84, 81, 59, 65, 47, 51, 157, 58];
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

%% Extract two-way inter-chromosome contacts
[~, ~, readIDReIdx] = unique(porecMatrix(:, 1)); % Reindex the ReadID
porecMatrix(:, 1) = readIDReIdx;
allChrContacts = cell(readIDReIdx(end), 1);
for i = 1:readIDReIdx(end)
    chrContact = unique(porecMatrix(porecMatrix(:, 1)==i, [2 4]));
    allChrContacts{i} = chrContact(:)';
end
twoWayInter = [20 22];
interChrLength = find((cellfun(@(v) all(ismember(v, twoWayInter)) & length(v)==length(twoWayInter), allChrContacts))>0);
interChrIntersect = intersect(porecMatrix(:, 1), interChrLength);
porecMatrixInterChr = porecMatrix(ismember(porecMatrix(:, 1), interChrIntersect), :);

%% Create an incidence matrix
[~, ~, readIDReIdx] = unique(porecMatrixInterChr(:, 1)); % Reindex the ReadID
porecMatrixInterChr(:, 1) = readIDReIdx;
incidenceMatrixA = sparse(porecMatrixInterChr(:, 3), porecMatrixInterChr(:, 1), 1,...
    sum(chrLength), porecMatrixInterChr(end, 1)); 
incidenceMatrixB = sparse(porecMatrixInterChr(:, 5), porecMatrixInterChr(:, 1), 1,...
    sum(chrLength), porecMatrixInterChr(end, 1));
incidenceMatrix = (incidenceMatrixA+incidenceMatrixB)>0;
incidenceMatrix = unique(incidenceMatrix', 'rows')'; % Remove duplicate edge

%% Extract all higer-order hyperedges
highOrderContacts = cell(size(incidenceMatrix, 2), 1);
for i = 1:size(incidenceMatrix, 2)
    highOrderContacts{i} = find(incidenceMatrix(:, i)>0);    
end
highOrderContacts = flipud(highOrderContacts(cellfun('length', highOrderContacts)>=3)); % Order >=3

%% Create CSV file for PAOHvis
CSVIncidenceMatrix = [];
label = [];
for i = 1:length(highOrderContacts)
    contact = highOrderContacts{i};
    for j = 1:length(contact)
        CSVIncidenceMatrix = [CSVIncidenceMatrix; [i contact(j)]]; %#ok<AGROW>
        if contact(j) >= 2725 && contact(j) <= 2789
            label{end+1} = 'Chr 20'; %#ok<*SAGROW>
        elseif contact(j) >= 2836 && contact(j) <= 2887
            label{end+1} = 'Chr 22';
        else 
            label{end+1} = 'Other';
        end
    end 
end 
        
CSVIncidenceMatrixTable = table(CSVIncidenceMatrix(:, 1), num2str(CSVIncidenceMatrix(:, 2), '%04d'), ...
    cell(size(CSVIncidenceMatrix, 1), 1), cell(size(CSVIncidenceMatrix, 1), 1), label');
writetable(CSVIncidenceMatrixTable, 'Inter_Chr_20_22.csv')

