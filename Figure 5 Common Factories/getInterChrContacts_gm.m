%% Load data
load('public_gm12878_reads.mat')
load('bad_locs_gm_100kb.mat')
binSize = 25e6;
% public_gm12878_strands = sortrows(public_gm12878_reads, 'read_id', 'ascend');
public_gm12878_reads.posA_25mb = ceil(public_gm12878_reads.posA/binSize);
public_gm12878_reads.posB_25mb = ceil(public_gm12878_reads.posB/binSize);
public_gm12878_reads.posA_100kb = ceil(public_gm12878_reads.posA/1e5);
public_gm12878_reads.posB_100kb = ceil(public_gm12878_reads.posB/1e5);
porecMatrix = table2array(public_gm12878_reads);

clear public_gm12878_reads

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
bad_read_locs = ismember(porecMatrix(:,8), find(bad_locs_gm_100kb))|...
                ismember(porecMatrix(:,9), find(bad_locs_gm_100kb));
            
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

%% Construct interChrContacts for most common (fig 5)
% output chr-chr
[~, ~, readIDReIdx] = unique(porecMatrix(:, 1)); % Reindex the ReadID
porecMatrix(:, 1) = readIDReIdx;
interChrSize = cell(readIDReIdx(end), 1);

parfor i = 1:readIDReIdx(end)
    x = unique(porecMatrix(porecMatrix(:, 1)==i, [2 4]));
    interChrSize{i} = x(:)';
end

interChrSizeChar = cellfun(@(v)sort(char(v)), interChrSize, 'uni', 0);
[interChrSizeUnique, ~, contactIdx] = unique(interChrSizeChar);
interChrSizeUnique = cellfun(@double, interChrSizeUnique, 'uni', 0);
weights = accumarray(contactIdx, 1);
interChrContacts = table(interChrSizeUnique, weights);    
    
save('gmInterChrContacts25mb.mat','interChrContacts')