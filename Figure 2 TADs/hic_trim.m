function [H_trim,bad_locs] = hic_trim(H,num_diag,num_sparse)
%hic_trim trims Hi-C data to remove poor mapping regions
%   H: Hi-C data matrix
%   num_diag: number of off diagonal elements to remove from consideration.
%
%   Default is to remove main diagonal from consideration ( =1)
%   num_sparse: number of bins in a column that must be non-zero. if >=1, 
%   this is the number of total bins in each column that must be >0.
%   if <1, this is the fraction of bins that must be >0.
%
%   Scott Ronquist, 2018

if nargin<2;num_diag=1;end
if nargin<3;num_sparse=0;end

% turn nan to zero
H(isnan(H)) = 0;

%% loop through samples and check for bad_locs
bad_locs = zeros(size(H,3),size(H,1));
for i = 1:size(H,3)
    % make sure symmetric
    H_temp = (H(:,:,i)+H(:,:,i)')./2;
    
    % remove diagonal
    H_temp = triu(H_temp-(triu(H_temp)-triu(H_temp,num_diag)));
    H_temp = H_temp+H_temp';
    
    % determine where matrix is zero
    bad_locs(i,:) = sum(H_temp)==0;
    
    % remove if sparse
    if num_sparse < 1
        temp = sum(logical(H_temp))./size(H_temp,1) < num_sparse;
        bad_locs(i,:) = logical(bad_locs(i,:)+temp);
    else
        temp = sum(logical(H_temp)) < num_sparse;
        bad_locs(i,:) = logical(bad_locs(i,:)+temp);
    end
    
end

%% trim matrix
%determine where there is a bad_loc in any sample and remove
bad_locs = logical(sum(bad_locs,1));

H_trim = H(~bad_locs,~bad_locs,:);

end

