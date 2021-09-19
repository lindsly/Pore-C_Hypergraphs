tic
% Add chromosome offsets
chrLength = [249, 243, 199, 191, 182, 171, 160, 146, 139, 134, 136,...
    134, 115, 108, 102, 91, 84, 81, 59, 65, 47, 51, 157, 58];

% If allChrContacts has been computed previously, load it and porecMatrix
if exist('allChrContacts_Fib.mat','file') == 2
    load('allChrContacts_Fib.mat')
else
    % Load dataset and convert to matrix
    load('v1234_np_uniqueIDs.mat')
    v1234_np_uniqueIDs = v1234_np_uniqueIDs(:, 2:end);
    binSize = 1e6;
    v1234_np_uniqueIDs.posA = ceil(v1234_np_uniqueIDs.posA/binSize);
    v1234_np_uniqueIDs.posB = ceil(v1234_np_uniqueIDs.posB/binSize);
    v1234_np_uniqueIDs = sortrows(v1234_np_uniqueIDs, 'read_id', 'ascend');
    porecMatrix = table2array(v1234_np_uniqueIDs);

    for chr = 2:24
        binAdd = sum(chrLength(1:chr-1));
        binIdxL = porecMatrix(:, 2)==chr;
        porecMatrix(binIdxL, 3) = porecMatrix(binIdxL, 3)+binAdd; 
        binIdxR = porecMatrix(:, 4)==chr;
        porecMatrix(binIdxR, 5) = porecMatrix(binIdxR, 5)+binAdd;
    end
    % Threshold contacts removing noises
    [~, uniqReadContacts, ~] = unique(porecMatrix(:, [1 3 5]), 'rows');
    porecMatrix = porecMatrix(uniqReadContacts, :);
    pairedContacts = porecMatrix(:, [3 5]);
    [~, ~, contactIdx] = unique(pairedContacts, 'rows', 'stable');  
    contactHist = accumarray(contactIdx, 1);
    contactCount = contactHist(contactIdx);
    eps = prctile(contactCount, 85);
    porecMatrix = porecMatrix(contactCount>=eps, :);

    %% Extract all inter-/intra- chromosome contacts
    [~, ~, readIDReIdx] = unique(porecMatrix(:, 1)); % Reindex the ReadID
    porecMatrix(:, 1) = readIDReIdx;
    allChrContacts = cell(readIDReIdx(end), 1);
    parfor i = 1:readIDReIdx(end)
        chrContact = unique(porecMatrix(porecMatrix(:, 1)==i, [2 4]));
        allChrContacts{i} = chrContact(:)';
    end
end

entropyChr = zeros(24, 24, 24);
entropyChrNorm = zeros(24, 24, 24);
nodeSize = zeros(24, 24, 24);
nodeContacts = zeros(24, 24, 24);
for chr1 = 1:24
    for chr2 = chr1:24
        for chr3 = chr2:24
            chrComb = unique([chr1, chr2, chr3]);
            interChrLength = find((cellfun(@(v) all(ismember(v, chrComb)) & length(v)==length(chrComb), allChrContacts))>0);
            interChrIntersect = intersect(porecMatrix(:, 1), interChrLength);
            porecMatrixInterChr = porecMatrix(ismember(porecMatrix(:, 1), interChrIntersect), :);
            
            if ~isempty(porecMatrixInterChr)
                % Create an incidence matrix
                [~, ~, readIDReIdx] = unique(porecMatrixInterChr(:, 1)); % Reindex the ReadID
                porecMatrixInterChr(:, 1) = readIDReIdx;
                incidenceMatrixA = sparse(porecMatrixInterChr(:, 3), porecMatrixInterChr(:, 1), 1,...
                    sum(chrLength), porecMatrixInterChr(end, 1));
                incidenceMatrixB = sparse(porecMatrixInterChr(:, 5), porecMatrixInterChr(:, 1), 1,...
                    sum(chrLength), porecMatrixInterChr(end, 1));
                incidenceMatrix = (incidenceMatrixA+incidenceMatrixB)>0;
                incidenceMatrix(all(~incidenceMatrix, 2), :) = [];
                incidenceMatrix = unique(incidenceMatrix', 'rows')'; % Remove duplicate edge

                % Construct Laplacian matrix            
                lapMatrix = incidenceMatrix*incidenceMatrix';
                eigValues = eig(lapMatrix);
                eigValues = eigValues(eigValues>0);
                eigValues = eigValues/sum(eigValues);
                entropyChr(chr1, chr2, chr3) = -sum(eigValues.*log(eigValues));
                nodeSize(chr1, chr2, chr3) = size(incidenceMatrix, 1);
                nodeContacts(chr1, chr2, chr3) = size(incidenceMatrix, 2);
                entropyChrNorm(chr1, chr2, chr3) = -sum(eigValues.*log(eigValues))/size(incidenceMatrix, 1);
            end
        end
    end
end
toc

