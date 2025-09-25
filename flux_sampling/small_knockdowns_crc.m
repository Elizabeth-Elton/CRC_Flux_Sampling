%% small script to run just the model knockdowns that I want in CRC model (adapted from Tavakoli et al 2025)

function [kd_model, minFlux, maxFlux] = small_knockdowns(j, percent, k, filename, save_model, fva, sample)
% gene_knockdowns.m runs the reaction for a given j (knockdown) at a given percent for a specific model in the model list and
% returns the knocked-down model

% input(s): j(int) - knockdown reaction number, percent(double) - percent
% to knockdown by (ie, percent - 0.2 means knocking down by 20%), k(int)
% model list index - this may change/not be needed

% output(s): kd_model(struct) - the knocked down modelt ouse in fva, flux
% sampling


% Load in the baseline data, calculated from original flux balance
% analysis. 
load('baselineData/KRAS_CRC_media_vs_CAF_media/data.mat');
KRAS_CRC=data(:,1:2:end);
%KRAS_CCM=data(:,2:2:end);

% Run the enzyme knockdowns. 

    disp(['Condition: ' 'KRAS CRC']);
    data_path = ['baselineData/KRAS_CRC_media_vs_CAF_media/data.mat'];
    load(data_path);

    b_nonzero_indices_mat = zeros(1, 1, 74, 89);  % 4 conditions, 100 sets of mass balance constraints, 74 reactions, 89 metabolites
    knockdown_biomass_mat = cell(1, 1, 1);

    model = WT_Model.model_lst{1, k}; % not actually the wt model, just the CRC is labelled as the wt model in this mat file
    model.b(82) = WT_Model.v_lst(73, k);

    model.b(72) = 0;

        %if statements that give fractional enzyme knockouts depending on
        %condition being used
        
            if mean(KRAS_CRC(:,j)) > 0
                    tmp_model = changeRxnBounds(model, model.rxns{j},(mean(KRAS_CRC(:,j)))*(1-percent), 'u');
                   % tmp_sol = optimizeCbModel(tmp_model, 'max');
                    tmp_sol = minimizeModelFlux(tmp_model, 'min', 'one');
            elseif  mean(KRAS_CRC(:,j)) < 0
                    tmp_model = changeRxnBounds(model, model.rxns{j},(mean(KRAS_CRC(:,j)))*(1-percent), 'l');
                    %tmp_sol = optimizeCbModel(tmp_model, 'max');
                    tmp_sol = minimizeModelFlux(tmp_model, 'min', 'one');
            elseif  mean(KRAS_CRC(:,j)) == 0
                    tmp_model = changeRxnBounds(model, model.rxns{j},(mean(KRAS_CRC(:,j)))*(1-percent), 'u');
                   % tmp_sol = optimizeCbModel(tmp_model, 'max');
                    tmp_sol = minimizeModelFlux(tmp_model, 'min', 'one');

            end

            if ~tmp_sol.stat
                relaxOption.internalRelax = 0;
                relaxOption.exchangeRelax = 0;
                relaxOption.steadyStateRelax = 1;
                relaxOption.epsilon = 10e-6;
                mets_to_exclude = false(size(tmp_model.mets));          
                mets_to_exclude(contains(tmp_model.mets, "_greater")) = true;
                mets_to_exclude(contains(tmp_model.mets, "_lower")) = true;
                relaxOption.excludedMetabolites = mets_to_exclude;
                solution = relaxedFBA(tmp_model,relaxOption);
                [stat,v,r,p,q] = deal(solution.stat, solution.v, solution.r, solution.p, solution.q);
                if stat == 0
                    v = zeros(size(tmp_model.rxns));
                    r = zeros(size(tmp_model.mets));
                    p = zeros(size(tmp_model.rxns));
                    q = zeros(size(tmp_model.rxns));
                    knockdown_biomass_mat{1, 1, j} = NaN;
                else
                    b_nonzero_indices = find(r ~= 0);
                    b_nonzero_values = r(b_nonzero_indices);
                    b_nonzero_indices_mat(1, 1, j, b_nonzero_indices) = b_nonzero_indices_mat(1, 1, j, b_nonzero_indices)+1;
                    for m = 1:numel(b_nonzero_indices)
                        tmp_model.b(b_nonzero_indices(m)) = tmp_model.b(b_nonzero_indices(m)) - b_nonzero_values(m);
                    end    
                   % solutionDel = optimizeCbModel(tmp_model, 'max');
                    solutionDel = minimizeModelFlux(tmp_model, 'min', 'one');
                    stat = solutionDel.stat;
                    if stat == 1
                       v = solutionDel.x;                      
                       if solutionDel.f > 0
                           knockdown_biomass_mat{1, 1, j} = solutionDel.v;
                       else
                           knockdown_biomass_mat{1, 1, j} = 0;
                       end
                    else
                       v = zeros(size(tmp_model.rxns));
                       r = zeros(size(tmp_model.mets));
                       p = zeros(size(tmp_model.rxns));
                       q = zeros(size(tmp_model.rxns));
                       knockdown_biomass_mat{1, 1, j} = 0;
                    end
                end
            else
                if tmp_sol.f > 0
                    knockdown_biomass_mat{1, 1, j} = tmp_sol.v;
                else
                    knockdown_biomass_mat{1, 1, j} = 0;
                end
            end

        if save_model == 0
            writeCbModel(tmp_model, 'fileName', strcat(filename, '.mat'))
        end
        
        tmp_modelz = changeCOBRAConstraints(tmp_model, 'Constraint1', 'dsense', 'G'); %SF says G is ok
        tmp_modelz.c = 0 * tmp_modelz.c; 
        
        if fva == 0
            [minFlux, maxFlux] = fluxVariability(tmp_modelz, [], 'min', [], 1, 1, 'heuristics', 0);
        end

        kd_model = tmp_model;

        if sample == 0
            options.nPointsReturned = 2000;
            options.nFiles = 1;
            [sampled_model, ~] =  sampleCbModel(tmp_modelz, strcat(filename, '_samples'), 'ACHR', options);

            writeCbModel(sampled_model, strcat(filename, '_achr.mat'))
        end
end

