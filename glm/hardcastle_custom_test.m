% Run this function with a list of cells, while in directory containing all
% day folders, to classify all cells in the given list. Outputs a txt file
% containing the classification of the given cells, as well as the p-values
% of all significance tests run for different models on the cells.

function hardcastle_custom_test(tbin_size, fc, filename, redo)
    % PARAMETERS:
    % tbin_size - size of time bin (in seconds) for binning of vmpv data.
    % fc - number of folds for cross-validation in glm_hardcastle.
    % cell_list - txt file containing list of cells to process.
    
    if ~exist('redo', 'var')
        redo = false;
    end

    % Read in list of cells from txt file
    if ~exist('filename', 'var')
        filename = 'cell_list.txt';  % should be in the directory level containing all day folders
    end
    cell_list = textread(filename, '%s', 'delimiter', '\n');
    
    curr_dir = pwd;
    for i = 1:length(cell_list)
        % Go into each cell directory, and run hardcastle functions
        cd([curr_dir '/' cell_list{i}]);
        disp(['Running classfication for ', pwd]);
        if exist('50ms_10fold_basetest(23-09-25)/glm_hardcastle_results.mat', 'file')
            load('50ms_10fold_basetest(23-09-25)/glm_hardcastle_results.mat', 'hc_results');
            hardcastle_cell_test(tbin_size, fc, redo, hc_results.params_consol);
        else
            hardcastle_cell_test(tbin_size, fc, redo);
        end
    end
    cd(curr_dir);
end
