%% Load data
% load('NlaIII_digest_IRFibroblast_006_reads.mat')
if ~exist('v1234_np_uniqueIDs','var')
    load('../IR_Fibroblast/v1234_np_uniqueIDs.mat')
    v1234_np_uniqueIDs(:,1) = [];
end
load('bad_locs_fib_100kb.mat')
binSize = 25e6;
v1234_np_uniqueIDs = sortrows(v1234_np_uniqueIDs, 'read_id', 'ascend');
v1234_np_uniqueIDs.posA_25mb = ceil(v1234_np_uniqueIDs.posA/binSize);
v1234_np_uniqueIDs.posB_25mb = ceil(v1234_np_uniqueIDs.posB/binSize);
v1234_np_uniqueIDs.posA_100kb = ceil(v1234_np_uniqueIDs.posA/1e5);
v1234_np_uniqueIDs.posB_100kb = ceil(v1234_np_uniqueIDs.posB/1e5);

porecMatrix = table2array(v1234_np_uniqueIDs);

% %%
% load('public_gm12878_reads.mat')
% binSize = 25e6;
% public_gm12878_reads.posA = ceil(public_gm12878_reads.posA/binSize);
% public_gm12878_reads.posB = ceil(public_gm12878_reads.posB/binSize);
% public_gm12878_strands = sortrows(public_gm12878_reads, 'read_id', 'ascend');
% porecMatrix = table2array(public_gm12878_reads);

%% Adjust bins
chrLength = ceil([249, 243, 199, 191, 182, 171, 160, 146, 139, 134, 136,...
    134, 115, 108, 102, 91, 84, 81, 59, 65, 47, 51, 157, 58]./25);
chrLength_100kb = ceil([248956422, 242193529, 198295559, 190214555, 181538259, 170805979, ...
    159345973, 145138636, 138394717, 133797422, 135086622, ...
    133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, ...
    58617616, 64444167, 46709983, 50818468, 156040895, 57227415]./1e5);
chrSum = zeros(1, 24);

for chr = 2:24
    binAdd = sum(chrLength(1:chr-1));
    binAdd_100kb = sum(chrLength_100kb(1:chr-1));
    binIdxL = porecMatrix(:, 2)==chr;
    porecMatrix(binIdxL, 6) = porecMatrix(binIdxL, 6)+binAdd; 
    porecMatrix(binIdxL, 8) = porecMatrix(binIdxL, 8)+binAdd_100kb; 
    binIdxR = porecMatrix(:, 4)==chr;
    porecMatrix(binIdxR, 7) = porecMatrix(binIdxR, 7)+binAdd;
    porecMatrix(binIdxR, 9) = porecMatrix(binIdxR, 9)+binAdd_100kb;

    
    chrSum(chr) = sum(chrLength(1:chr-1));
end
chrSum = [chrSum, 136];

%% Remove outlier reads near centromeres
bad_read_locs = ismember(porecMatrix(:,8), find(bad_locs_fib_100kb))|...
                ismember(porecMatrix(:,9), find(bad_locs_fib_100kb));
            
porecMatrix(bad_read_locs,:) = [];


%% Threshold contacts removing noise
[~, uniqReadContacts, ~] = unique(porecMatrix(:, [1 6 7]), 'rows');
porecMatrix = porecMatrix(uniqReadContacts, :);
pairedContacts = porecMatrix(:, [6 7]);
[~, ~, contactIdx] = unique(pairedContacts, 'rows', 'stable');  
contactHist = accumarray(contactIdx, 1);
contactCount = contactHist(contactIdx);
eps = prctile(contactCount, 85);
porecMatrix = porecMatrix(contactCount>=eps, :);

%% Create an incidence matrix
[~, ~, readIDReIdx] = unique(porecMatrix(:, 1)); % Reindex the ReadID
porecMatrix(:, 1) = readIDReIdx;

incidenceMatrixA = sparse(porecMatrix(:, 6), porecMatrix(:, 1), 1,...
    sum(chrLength), porecMatrix(end, 1)); 
incidenceMatrixB = sparse(porecMatrix(:, 7), porecMatrix(:, 1), 1,...
    sum(chrLength), porecMatrix(end, 1));
incidenceMatrix = (incidenceMatrixA+incidenceMatrixB)>0;
incidenceMatrix = incidenceMatrix(:, sum(incidenceMatrix, 1)>2);
[incidenceMatrix, ~, contactIdx] = unique(incidenceMatrix', 'rows', 'stable');
incidenceMatrix = incidenceMatrix';

%% Debugging centromeric reads
% Row sum of incidence matrix for degree per locus
degree_cent_fib_25mb = sum(incidenceMatrix,2);
% bad_locs_fib_100kb = degree_cent_fib_100kb>1000;

% 
% figure, hist(degree_cent_fib_100kb,1000)
% prctile(degree_cent_fib_100kb,99.95)
% toc


% Extract incidence Matrix with Chromosome i
weights = accumarray(contactIdx, 1);
% Chr 24 (Y) is removed during filtering noise step
incidenceMatrixChosen = zeros(size(incidenceMatrix, 1), 10*23);
for i = 1:23
    incidenceMatrixChr = incidenceMatrix(chrSum(i)+1:chrSum(i+1), :);
    weightsChr = weights(sum(incidenceMatrixChr, 1)>0);
    incidenceMatrixChr = incidenceMatrix(:, sum(incidenceMatrixChr, 1)>0);         
    
    
    incidenceMatrixChrIntra = incidenceMatrixChr(:, sum(incidenceMatrixChr, 1)==sum(incidenceMatrixChr(chrSum(i)+1:chrSum(i+1), :), 1));
    weightsIntra = weightsChr(sum(incidenceMatrixChr, 1)==sum(incidenceMatrixChr(chrSum(i)+1:chrSum(i+1), :), 1));
    
    incidenceMatrixChrInter = incidenceMatrixChr(:, ~ismember(incidenceMatrixChr', incidenceMatrixChrIntra', 'rows'));
    weightsInter = weightsChr(~ismember(incidenceMatrixChr', incidenceMatrixChrIntra', 'rows'));
    
    incidenceMatrixChrTemp = zeros(size(incidenceMatrix, 1), 10);
    if size(incidenceMatrixChrIntra, 2) >= 5
        [~, weightsIntraIdx] = sort(weightsIntra, 'descend');
        [~, weightsInterIdx] = sort(weightsInter, 'descend');
        incidenceMatrixChrTemp(:, 1:5) = incidenceMatrixChrIntra(:, weightsIntraIdx(1:5));
        incidenceMatrixChrTemp(:, 6:10) = incidenceMatrixChrInter(:, weightsInterIdx(1:5));
    elseif size(incidenceMatrixChrIntra, 2) < 5 && ~isempty(incidenceMatrixChrIntra)
        [~, weightsIntraIdx] = sort(weightsIntra, 'descend');
        [~, weightsInterIdx] = sort(weightsInter, 'descend');
        incidenceMatrixChrTemp(:, 1:size(incidenceMatrixChrIntra, 2)) = incidenceMatrixChrIntra(:, weightsIntraIdx(1:size(incidenceMatrixChrIntra, 2)));
        incidenceMatrixChrTemp(:, size(incidenceMatrixChrIntra, 2)+1:10) = incidenceMatrixChrInter(:, weightsInterIdx(1:10-size(incidenceMatrixChrIntra, 2)));
    else
        [~, weightsInterIdx] = sort(weightsInter, 'descend');
        incidenceMatrixChrTemp = incidenceMatrixChrInter(:, weightsInterIdx(1:10));
    end 
    incidenceMatrixChosen(:, (i-1)*10+1:i*10) = incidenceMatrixChrTemp;
    weights(ismember(incidenceMatrix', incidenceMatrixChrTemp', 'rows')) = [];
    incidenceMatrix(:, ismember(incidenceMatrix', incidenceMatrixChrTemp', 'rows')) = [];
end

%% Output edges
edgeSet = cell(size(incidenceMatrixChosen, 2), 1);
for i = 1:size(incidenceMatrixChosen, 2)
    edgeSet{i} = find(incidenceMatrixChosen(:, i)>0)';
end

degree_chosen_fib_25mb = sum(incidenceMatrixChosen,2);


%% Create .csv file for Paovis
CSVIncidenceMatrix = [];
label = [];
for i = 1:length(edgeSet)
    contact = edgeSet{i};
    for j = 1:length(contact)
        CSVIncidenceMatrix = [CSVIncidenceMatrix; [i contact(j)]]; %#ok<AGROW>
        chrDiff = contact(j)-chrSum;
        chrDiff = chrDiff(chrDiff>0);
        label = [label; length(chrDiff)]; %#ok<AGROW>
    end
end

CSVIncidenceMatrixTable = table(CSVIncidenceMatrix(:, 1), num2str(CSVIncidenceMatrix(:, 2), '%03d'), ...
    cell(size(CSVIncidenceMatrix, 1), 1), cell(size(CSVIncidenceMatrix, 1), 1), label);

writetable(CSVIncidenceMatrixTable, 'CSVIncidenceMatrixTableFib.csv')    


    
    
