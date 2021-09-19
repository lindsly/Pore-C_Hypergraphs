%% Load data
type = 2;
if type == 1
    load('NlaIII_digest_IRFibroblast_006_reads.mat')
    binSize = 1e5;
    NlaIII_digest_IRFibrobalst_006_reads = sortrows(NlaIII_digest_IRFibrobalst_006_reads, 'read_id', 'ascend');
    NlaIII_digest_IRFibrobalst_006_reads.posA = ceil(NlaIII_digest_IRFibrobalst_006_reads.posA/binSize);
    NlaIII_digest_IRFibrobalst_006_reads.posB = ceil(NlaIII_digest_IRFibrobalst_006_reads.posB/binSize);
    porecMatrix = table2array(NlaIII_digest_IRFibrobalst_006_reads);
elseif type == 2
    load('public_gm12878_reads.mat')
    binSize = 1e5;
    public_gm12878_strands = sortrows(public_gm12878_reads, 'read_id', 'ascend');
    public_gm12878_reads.posA = ceil(public_gm12878_reads.posA/binSize);
    public_gm12878_reads.posB = ceil(public_gm12878_reads.posB/binSize);
    porecMatrix = table2array(public_gm12878_reads);
end
%% Adjust bins
chrLength = ceil([248956422, 242193529, 198295559, 190214555, 181538259, 170805979, ...
    159345973, 145138636, 138394717, 133797422, 135086622, ...
    133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, ...
    58617616, 64444167, 46709983, 50818468, 156040895, 57227415]./binSize);

for chr = 2:24
    binAdd_100kb = sum(chrLength(1:chr-1));
    binIdxL = porecMatrix(:, 2)==chr;
    porecMatrix(binIdxL, 3) = porecMatrix(binIdxL, 3)+binAdd_100kb; 
    binIdxR = porecMatrix(:, 4)==chr;
    porecMatrix(binIdxR, 5) = porecMatrix(binIdxR, 5)+binAdd_100kb;
end

%% Threshold contacts removing noise
[~, uniqReadContacts, ~] = unique(porecMatrix(:, [1 3 5]), 'rows');
porecMatrix = porecMatrix(uniqReadContacts, :);
pairedContacts = porecMatrix(:, [3 5]);
[~, ~, contactIdx] = unique(pairedContacts, 'rows', 'stable');  
contactHist = accumarray(contactIdx, 1);
contactCount = contactHist(contactIdx);
eps = prctile(contactCount, 85);
porecMatrix = porecMatrix(contactCount>=eps, :);

%% Create an incidence matrix
[~, ~, readIDReIdx] = unique(porecMatrix(:, 1)); % Reindex the ReadID
porecMatrix(:, 1) = readIDReIdx;

incidenceMatrixA = sparse(porecMatrix(:, 3), porecMatrix(:, 1), 1,...
    sum(chrLength), porecMatrix(end, 1)); 
incidenceMatrixB = sparse(porecMatrix(:, 5), porecMatrix(:, 1), 1,...
    sum(chrLength), porecMatrix(end, 1));
incidenceMatrix = (incidenceMatrixA+incidenceMatrixB)>0;
incidenceMatrix = incidenceMatrix(:, sum(incidenceMatrix, 1)>2);
incidenceMatrix = unique(incidenceMatrix', 'rows')';

%% Debugging centromeric reads
% Row sum of incidence matrix for degree per locus to find distribution
if type == 1
    degree_cent_fib_100kb = sum(incidenceMatrix,2);
    figure, histogram(degree_cent_fib_100kb,1000)
    bad_locs_fib_100kb = degree_cent_fib_100kb>prctile(degree_cent_fib_100kb,99.9);
%     prctile(degree_cent_fib_100kb,99.9)
elseif type == 2
    degree_cent_gm_100kb = sum(incidenceMatrix,2);
    figure, histogram(degree_cent_gm_100kb,1000)
    bad_locs_gm_100kb = degree_cent_gm_100kb > prctile(degree_cent_gm_100kb,99.9);
end

% 
% figure, hist(degree_cent_fib_100kb,1000)
% prctile(degree_cent_fib_100kb,99.95)
% toc

% bad_read_locs