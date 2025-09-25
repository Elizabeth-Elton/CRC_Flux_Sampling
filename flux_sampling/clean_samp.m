%  perform fva and ACHR sampling on base model and knockdowns FOR CCM
%  CONDITION

%addpath('/Users/.../cobratoolbox')
initCobraToolbox(false);
%%
changeCobraSolver('gurobi');
%% get the ccm the base model 
load('flux_sampling/baselineData/KRAS_CRC_media_vs_CAF_media/data.mat');

model = KRAS_Model.model_lst{1, 1};

modelz = changeCOBRAConstraints(model, 'Constraint1', 'dsense', 'G'); %G is ok - not strictly necessary but very few things (if any) vary if i do not

modelz.c = 0 * modelz.c; %clear objective function 

%fva
[minFs3, maxFs3] = fluxVariability(modelz, [], 'min', [], 2, 1, 'heuristics', 0);
%% sample base model
options.nPointsReturned = 2000;
options.nFiles = 1;
[base_model, base_samples] =  sampleCbModel(modelz, [], 'ACHR', options);

%save the base model to use later
writeCbModel(base_model, 'base_model_achr.mat')

%% get list of top 5% of models

opts = detectImportOptions("data/top_kd_rxns_ccm.csv");
opts.Delimiter = ',';
opts.VariableNamingRule = 'preserve';
opts.VariableNames = {'rxn_number', 'percent_kd', 'reaction_name'};
opts = setvartype(opts,{'rxn_number', 'percent_kd'},'double');
opts = setvartype(opts,{'reaction_name'},'char');

rxn_inputs = readtable('data/top_kd_rxns_ccm.csv', opts);
%% sample the top models

modelz = cell(1, 12);
minFluxez = cell(1, 12);
maxFluxez = cell(1, 12);

for i=1:12
    [mod1, minf1, maxf1] = small_knockdowns(rxn_inputs.rxn_number(i), rxn_inputs.percent_kd(i),... 
        1, rxn_inputs.reaction_name{i}, 0, 0, 0);
    modelz{i} = mod1;
    minFluxez{i} = minf1;
    maxFluxez{i} = maxf1;

end

%% compare the reaction lists for each of these

% get all of those files i just created
folder_path = 'data/sampling_results';

% separate out into models and sampling results
model_list = dir(fullfile(folder_path, ['*' 'achr' '*.mat']));
samples_list = dir(fullfile(folder_path, ['*' 'samples' '*.mat']));

all_my_models = cell(1, length(model_list));
all_my_samples = cell(1, length(samples_list));

% loop through all the model files
for i = 1:length(model_list)

    file_path = fullfile(folder_path, model_list(i).name);
    loaded_data = load(file_path);
    all_my_models{i} = loaded_data;
   
end

% loop through all the sampling files 
for i = 1:length(samples_list)

    file_path = fullfile(folder_path, samples_list(i).name);
    loaded_data = load(file_path);
    all_my_samples{i} = loaded_data;
   
    fprintf('Loaded file: %s\n', samples_list(i).name);
end

%% match the reaction names with the sampling rows

labelled_samples = cell(1, 12);

for i = 1:12

    samples = all_my_samples{i}.points;
    r_names = all_my_models{i}.model.rxnNames;
    nice_data = horzcat(r_names, num2cell(samples));

    labelled_samples{i} = nice_data;

end
%% load in base model data and prep for kl divergence calc

base_points = load("_1.mat"); %base_samples.points
base_model = load('base_model_achr.mat');

base_model_names = base_model.model.rxnNames;
base_model_points = base_points.points;

nice_base_points = horzcat(base_model_names, num2cell(base_model_points));

divergences = cell(1,12);
labels = cell(1,12);
no_matches = zeros(1, 12);

%%
for i = 1:length(labelled_samples) %length of labelled samples 

    base_samp_to_edit = nice_base_points;
    kd_samp_i_want = labelled_samples{1,i};

    idx = find(ismember(base_samp_to_edit(:,1), kd_samp_i_want(:,1)) == 0); 

    no_match_count = length(idx); %total rxns that dont exist in the base model
    no_matches(i) = no_match_count;

    % where idx = 0, remove that rxn from base
    base_samp_to_edit(idx,:) = [];

    % check here if they dont match lengths
    if size(base_samp_to_edit, 1) ~= size(kd_samp_i_want, 1)
        print('Number of reactions does not match. check knockdown model reactions')
    end

    base_samp_to_edit(:,1) = []; % removing the labelled column so i can do cell2mat
    kd_samp_i_want(:,1) = [];

    %then do kl divergence
    distances = zeros(size(base_samp_to_edit, 1), 1);
    for j=1:size(base_samp_to_edit, 1) %loop through each reaction

        base_samp_edited = cell2mat(base_samp_to_edit);
        kd_samp_edited = cell2mat(kd_samp_i_want);

        [f1, x1] = ksdensity(base_samp_edited(j, :)); 
        [f2, x2] = ksdensity(kd_samp_edited(j, :));

        dist1 = KLDis(f1, f2);
        dist2 = KLDis(f2, f1);

        avg_dist = (dist1 + dist2)/2;

        distances(j) = avg_dist;

    end

    divergences{1,i} = distances;
    labels{1,i} = labelled_samples{1,i}(:,1);
end
%% export everything and move to r
teststruct = struct('divergence_vals', divergences, 'reaction_labels', labels);

save('data/model_divergences.mat', 'divergences')
save('data/sampling_results.mat', 'labelled_samples')
save('data/model_names.mat', 'model_list')

save('data/rxn_names_labelled.mat', 'labels') 

save('data/divergences_labelled', 'teststruct')