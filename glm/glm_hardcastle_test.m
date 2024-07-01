function hc_results = glm_hardcastle(glm_data, fc, smooth_params, init_params)
%
%   Reference: Hardcastle et al., 2017. A Multiplexed, Heterogeneous, and Adaptive Code for Navigation in Medial Entorhinal Cortex
%   Some code adapted from github.com/GiocomoLab/ln-model-of-mec-neurons
%
%   Be in directory of spiketrain.mat
%   Afterwards, do signrank tests, etc. with hardcastle_testing function
%   
%   PARAMETERS:
%   glm_data - uses either glm_vmpvData or glm_genData.
%   e.g. glm_hardcastle(glm_vmpvData(0.020), 10)
%
%   fc - number of folds used for cross-validation.
%
%   smooth_params - 1x3 array of beta values for smoothing parameters,
%   in the order of [ place, headdirection, view ]. Optional argument,
%   defaults to [ 3, 3, 3 ] if not provided.
%
%   init_params - fc x 7 array of pre-known optimized parameters to be used
%   for initialization, usually taken from a previous run on the same data

bin_stc = glm_data.bin_stc;
tbin_size = glm_data.tbin_size;
place_good_bins = glm_data.place_good_bins;
view_good_bins = glm_data.view_good_bins;
good_bins = { place_good_bins, view_good_bins };

% Define bin geometry of the environment, in terms of:
% 1. no of floor_width bins, 2. no of wall_height bins, 3. no of
% pillar_height bins, 4. no of hd bins
floor_width = 40; wall_height = 8; pillar_height = 5; hd_bins = 60;
bin_geom = [ floor_width, wall_height, pillar_height, hd_bins ];
% additional intermediate variables for viewspace bins
viewbin_offset = 2; wall_perim = floor_width*4; pillar_width = floor_width/5; pillar_perim = pillar_width*4;

first_feature_bins = floor_width^2; % default 1600
second_feature_bins = hd_bins; % default 60
third_feature_bins = viewbin_offset + 2*floor_width^2 + wall_height*wall_perim + 4*pillar_height*pillar_perim; % default 5122

% Beta values for smoothing parameters
if exist('smooth_params', 'var') && ~isempty(smooth_params)
    if length(smooth_params) == 1
        betas = [ smooth_params, smooth_params, smooth_params ]; % use same value for all variables if only 1 value is specified
    else
        betas = smooth_params;
    end
else
    betas = [ 3e0, 3e0, 3e0 ]; % [ beta_place, beta_hd, beta_view ], default value of 3 for each param
end

% Num shuffles and shuffle timestep limits for generation of null distribution
numShuffles = 50;
shuffleLimits = [0.1, 0.9];

% Start to fill in x and y matrices
samples_total = size(bin_stc,1);
x = zeros(samples_total,first_feature_bins+second_feature_bins+third_feature_bins);
y = bin_stc(1:end,5);

for k = 1:samples_total
    x(k,bin_stc(k,2)) = 1;
    x(k,first_feature_bins+bin_stc(k,3)) = 1;
    x(k,first_feature_bins+second_feature_bins+bin_stc(k,4)) = 1;
end

% Filters for unoccupied place and view bins
place_filter = ones(first_feature_bins,1);
view_filter = ones(third_feature_bins,1);
place_filter(place_good_bins) = 0;
view_filter(view_good_bins) = 0;
large_negative_number = -1e3;

%%% inputs done %%%

%%%%%%%%%% looping train-test splits %%%%%%%%%%

modelType = [1 1 1;
    1 1 0;
    1 0 1;
    0 1 1;
    1 0 0;
    0 1 0;
    0 0 1]; % testing for mixed model, ph model, pv model, hv model, place model, hd model, view model
modelName = {'phv', 'ph', 'pv', 'hv', 'place', 'headdirection', 'spatialview'};

folds = fc;
num_models = size(modelType, 1);

edges = round(linspace(1,length(y)+1, (5*folds)+1)); % splits dataset into 5xfold sections, to sample folds across time

testFit = nan(folds, num_models);
trainFit = nan(folds, num_models);
testFit_pure = nan(folds, num_models);
trainFit_pure = nan(folds, num_models);
paramsAll = cell(folds, num_models);
testSig = nan(folds, num_models);
trainSig = nan(folds, num_models);
testNullDist = cell(folds, num_models);
trainNullDist = cell(folds, num_models);

for model_type = 1:num_models % test different models on this dataset

    % Set bins that have no occurences to a large negative number during initialization of params
    bin_filter = [];
    if (modelType(model_type,1))
        bin_filter = [bin_filter; place_filter];
    end
    if (modelType(model_type,2))
        bin_filter = [bin_filter; ones(second_feature_bins,1)];
    end
    if (modelType(model_type,3))
        bin_filter = [bin_filter; view_filter];
    end
    bin_filter = find(bin_filter);
    
    if ~exist('init_params', 'var')
        % Random initialization of params for the first fold, then reuse
        % optimized params from the previous fold for subsequent folds
        param = 1e-3*randn(first_feature_bins*modelType(model_type,1) + second_feature_bins*modelType(model_type,2) + third_feature_bins*modelType(model_type,3), 1); % random initialization
        param(bin_filter) = large_negative_number; % set all bins that have no observations to -1e1 (sufficiently large negative number)
    end
    
    disp(['Fitting model ', modelName{model_type}])

    for k = 1:folds
        fprintf('Fold %d of %d \n', k, folds)
%       disp('selectivity, params used, fold');
%       disp([cell_type model_type k]);

        if exist('init_params', 'var')
            param = init_params{k, model_type};
        end

        test_ind  = [edges(k):edges(k+1)-1 edges(k+folds):edges(k+folds+1)-1 ...
            edges(k+2*folds):edges(k+2*folds+1)-1 edges(k+3*folds):edges(k+3*folds+1)-1 ...
            edges(k+4*folds):edges(k+4*folds+1)-1]; % grab indices for this fold, basically 5 smaller, spaced out sections

        train_ind = setdiff(1:length(y), test_ind); % remaining datapoints

        train_spikes = y(train_ind);
        test_spikes = y(test_ind);

        switch model_type
            case 1 % keep all columns (place, head direction and view info)
                train_A = x(train_ind,:);
                test_A = x(test_ind,:);
            case 2 % drop view info only
                train_A = x(train_ind,1:first_feature_bins+second_feature_bins);
                test_A = x(test_ind,1:first_feature_bins+second_feature_bins);
            case 3 % drop head direction info only
                train_A = x(train_ind,[1:first_feature_bins 1+first_feature_bins+second_feature_bins:end]);
                test_A = x(test_ind,[1:first_feature_bins 1+first_feature_bins+second_feature_bins:end]);
            case 4 % drop place info only
                train_A = x(train_ind,1+first_feature_bins:end);
                test_A = x(test_ind,1+first_feature_bins:end);
            case 5 % keep place info only
                train_A = x(train_ind,1:first_feature_bins);
                test_A = x(test_ind,1:first_feature_bins);
            case 6 % keep head direction info only
                train_A = x(train_ind,1+first_feature_bins:first_feature_bins+second_feature_bins);
                test_A = x(test_ind,1+first_feature_bins:first_feature_bins+second_feature_bins);
            case 7 % keep view info only
                train_A = x(train_ind,1+first_feature_bins+second_feature_bins:end);
                test_A = x(test_ind,1+first_feature_bins+second_feature_bins:end);            
        end

        opts = optimset('Gradobj','on','Hessian','on','Display','off');
        data{1} = train_A; 
        data{2} = train_spikes;
        init_param = param;
        % bottom part all adapted from reference github code
        [param] = fminunc(@(param) ln_poisson_model_vmpv(param,data,modelType(model_type,:),bin_geom,betas,good_bins), init_param, opts);

        % test fit

        r = exp(test_A * param); 
        n = test_spikes; 
        log_llh_test_model = nansum(r-n.*log(r)+log(factorial(n)))/sum(n); %note: log(gamma(n+1)) will be unstable if n is large (which it isn't here)
        
        % generate null distribution from shuffled timeseries
        log_llh_test_null = nan(numShuffles, 1);
        timeShifts = round(length(r) * (rand(numShuffles, 1) * diff(shuffleLimits) + shuffleLimits(1)));
        for i = 1:numShuffles
            r_shifted = circshift(r, timeShifts(i));
            log_llh_test_null(i) = nansum(r_shifted-n.*log(r_shifted)+log(factorial(n)))/sum(n);
        end
        
        testFit(k, model_type) = log(2)*(-log_llh_test_model + nanmean(log_llh_test_null));
        testFit_pure(k, model_type) = log(2)*(-log_llh_test_model);
        testSig(k, model_type) = signrank(log_llh_test_null, log_llh_test_model, 'tail', 'left');
        testNullDist{k, model_type} = log(2)*(-log_llh_test_null);

        % train fit

        r_train = exp(train_A * param); 
        n_train = train_spikes;
        log_llh_train_model = nansum(r_train-n_train.*log(r_train)+log(factorial(n_train)))/sum(n_train);
        
        % generate null distribution from shuffled timeseries
        log_llh_train_null = nan(numShuffles, 1);
        timeShifts = round(length(r_train) * (rand(numShuffles, 1) * diff(shuffleLimits) + shuffleLimits(1)));
        for i = 1:numShuffles
            r_shifted = circshift(r_train, timeShifts(i));
            log_llh_train_null(i) = nansum(r_shifted-n_train.*log(r_shifted)+log(factorial(n_train)))/sum(n_train);
        end
        
        trainFit(k, model_type) = log(2)*(-log_llh_train_model + nansum(log_llh_train_null));
        trainFit_pure(k, model_type) = log(2)*(-log_llh_train_model);
        trainSig(k, model_type) = signrank(log_llh_train_null, log_llh_train_model, 'tail', 'left');
        trainNullDist{k, model_type} = log(2)*(-log_llh_train_null);

        paramsAll{k, model_type} = param;

    end
    
end
    
%glm_hardcastle_results.inputs = [x y];
hc_results.training_fits = trainFit; % for each dataset (phv, ph, pv, hv, place, hd, view behavior), store n_folds x model_type training likelihood values (with subtraction of mean model)
hc_results.training_fits_pure = trainFit_pure; % training likelihood without comparison to null distribution
hc_results.training_sigs = trainSig; % significance testing of training likelihood against null distribution
hc_results.training_null_dist = trainNullDist; % null distribution of training log llh values
hc_results.testing_fits = testFit; % test likelihood with comparison to null distribution
hc_results.testing_fits_pure = testFit_pure; % test likelihood without comparison to null distribution
hc_results.testing_sigs = testSig; % significance testing of test likelihood against null distribution
hc_results.testing_null_dist = testNullDist; % null distribution of testing log llh values
hc_results.params_consol = paramsAll; % weights stored here

hc_results.tbin_size = tbin_size;
hc_results.num_folds = fc;
hc_results.smoothing_beta = betas;
%hc_results.ThresVel = glm_data.ThresVel;
%hc_results.UseMinObs = glm_data.UseMinObs;

save('glm_hardcastle_results.mat','hc_results','-v7.3');

end
