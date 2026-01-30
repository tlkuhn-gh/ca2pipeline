classdef NetworkBursts
    methods(Static)
        function network_bursts = detect_network_bursts(calcium_event_times, min_neurons_required, burst_frame_range)
        % DETECT_NETWORK_BURSTS
        % Identifies network bursts based on simultaneous calcium events.
        % 
        % Based on ISIn method from Bakkum, D. J., Radivojevic, M., Frey, U., 
        %   Franke, F., Hierlemann, A., & Takahashi, H. (2014). 
        %   Parameters for burst detection. Frontiers in computational neuroscience, 7, 193. 
        %   https://doi.org/10.3389/fncom.2013.00193
        %
        % Inputs:
        %   calcium_event_times : Nx1 cell array
        %       Each cell contains a vector of frame indices where that neuron fired.
        %   min_neurons_required: integer
        %       Minimum number of neurons required to define a network burst.
        %   burst_frame_range : integer
        %       Maximum frame difference allowed within a single burst.
        %       NOTE: This corresponds to the 'Synchronized time window' parameter
        %       in MEAToolbox's Network Bursts feature. MEAToolbox's default value is .1 s, 
        %       equivalent to 2 frames at 30 Hz.  
        %
        % Output:
        %   network_bursts : vector of doubles
        %       First frame index of each detected network burst.
        
            num_neurons = numel(calcium_event_times);
        
            % Collect all events with neuron identity
            all_frames = [];
            all_neuron_ids = [];
        
            for n = 1:num_neurons
                frames = calcium_event_times{n}(:);
                all_frames = [all_frames; frames];
                all_neuron_ids = [all_neuron_ids; repmat(n, numel(frames), 1)];
            end
        
            if isempty(all_frames)
                network_bursts = [];
                return;
            end
        
            % Sort events by frame
            [all_frames, sort_idx] = sort(all_frames);
            all_neuron_ids = all_neuron_ids(sort_idx);
        
            network_bursts = [];
            burst_start_idx = 1;
        
            while burst_start_idx <= numel(all_frames)
                burst_start_frame = all_frames(burst_start_idx);
        
                % Find all events within burst_frame_range
                burst_end_idx = burst_start_idx;
                while burst_end_idx < numel(all_frames) && ...
                      all_frames(burst_end_idx + 1) - burst_start_frame <= burst_frame_range
                    burst_end_idx = burst_end_idx + 1;
                end
                        
                % Use this line if using ISIn approach from Bakkum, et al.,
                % 2014 (need N spikes within ISIn ms, spikes from same neuron count towards N, N should be integer not percent of neurons)
                participating_neurons = all_neuron_ids(burst_start_idx:burst_end_idx);
        
                if numel(participating_neurons) >= min_neurons_required
                    network_bursts(end+1,1) = burst_start_frame; %#ok<AGROW>

                    % Check if following frames should be included in burst
                    % "The burst ends when this condition [N spikes in less
                    % than T ms] is no longer met." Bakkum, et al., 2014
                    while numel(participating_neurons) >= min_neurons_required
                        burst_start_idx = burst_start_idx + 1
                        burst_start_frame = all_frames(burst_start_idx);

                        burst_end_idx = burst_start_idx;
                        while burst_end_idx < numel(all_frames) && ...
                              all_frames(burst_end_idx + 1) - burst_start_frame <= burst_frame_range
                            burst_end_idx = burst_end_idx + 1;
                        end

                        participating_neurons = all_neuron_ids(burst_start_idx:burst_end_idx);
                    end

                    % Move to next candidate burst
                    burst_start_idx = burst_end_idx;
                else
                    burst_start_idx = burst_start_idx + 1;
                end
            end
        end
        
        function plot_network_burst_raster(calcium_event_times, network_bursts, graph_title)
        % PLOT_NETWORK_BURST_RASTER
        % Raster plot of calcium events with network burst markers.
        
            if nargin < 3
                graph_title = '';
            end
        
            num_neurons = numel(calcium_event_times);
        
            % Determine maximum frame index
            max_frame = 0;
            for i = 1:num_neurons
                if ~isempty(calcium_event_times{i})
                    max_frame = max(max_frame, max(calcium_event_times{i}));
                end
            end
        
            figure;
            ax = axes;
            hold(ax, 'on');
        
            % --- Plot network bursts FIRST (background) ---
            for i = 1:numel(network_bursts)
                xline(ax, network_bursts(i), ...
                    'r-', 'LineWidth', 2);
            end
        
            % --- Plot calcium events SECOND (foreground) ---
            for neuron_idx = 1:num_neurons
                frames = calcium_event_times{neuron_idx};
                if ~isempty(frames)
                    y_vals = repmat(neuron_idx, size(frames));
                    plot(ax, frames, y_vals, '.', ...
                        'Color', 'k', ...
                        'MarkerSize', 10);
                end
            end
        
            % Axis formatting
            xlabel(ax, 'Frame');
            ylabel(ax, 'Neuron Index');
            title(ax, graph_title);
        
            ylim(ax, [0.5, num_neurons + 0.5]);
            xlim(ax, [0, max_frame]);
            set(ax, 'YDir', 'reverse');  % Neuron 1 at top
            box(ax, 'on');
        
            hold(ax, 'off');
        end
        
        
        function batch_network_bursts = batch_process_network_bursts( ...
            root_folder, include_subfolders, file_prefix, ...
            percent_neurons, burst_frame_range, output_file_name)
        % BATCH_PROCESS_NETWORK_BURSTS
        % Batch-detects network bursts across multiple .mat files.
        %
        % Inputs:
        %   root_folder : string or char
        %       Folder to search for .mat files
        %   include_subfolders : logical (default = false)
        %       Whether to recursively search subfolders
        %   file_prefix : string or char
        %       Only process .mat files beginning with this prefix
        %   percent_neurons : scalar (0–100)
        %       Minimum percentage of neurons for a network burst
        %   burst_frame_range : integer
        %       Maximum frame difference for burst grouping
        %   output_file_name : string or char
        %       Base name (no extension) for output files
        %
        % Output:
        %   batch_network_bursts : table
        %       Summary table of detected network bursts
        
            % Default for include_subfolders
            if nargin < 2 || isempty(include_subfolders)
                include_subfolders = false;
            end
        
            % Build file search pattern
            if include_subfolders
                file_pattern = fullfile(root_folder, '**', [file_prefix '*.mat']);
            else
                file_pattern = fullfile(root_folder, [file_prefix '*.mat']);
            end
        
            files = dir(file_pattern);
        
            % Preallocate table variables
            FileName = {};
            FullFilePath = {};
            TotalNetworkBursts = [];
            BurstFrameIndices = {};
        
            for i = 1:numel(files)
                file_name = files(i).name;
                full_path = fullfile(files(i).folder, file_name);
        
                try
                    S = load(full_path);
                catch
                    fprintf('ERROR loading file: %s\n', full_path);
                    continue;
                end
        
                % Validate required structure
                if ~isfield(S, 'data') || ~isfield(S.data, 'Spikes_cell')
                    fprintf('ERROR: File %s does not contain data.Spikes_cell\n', full_path);
                    continue;
                end
        
                % Detect network bursts
                network_bursts = NetworkBursts.detect_network_bursts( ...
                    S.data.Spikes_cell, percent_neurons, burst_frame_range);
        
                % Strip prefix and .mat from filename
                stripped_name = erase(file_name, file_prefix);
                stripped_name = erase(stripped_name, '.mat');
        
                % Store results
                FileName{end+1,1} = stripped_name; %#ok<AGROW>
                FullFilePath{end+1,1} = full_path; %#ok<AGROW>
                TotalNetworkBursts(end+1,1) = numel(network_bursts); %#ok<AGROW>
                BurstFrameIndices{end+1,1} = network_bursts; %#ok<AGROW>
            end
        
            % Create output table
            batch_network_bursts = table( ...
                FileName, FullFilePath, TotalNetworkBursts, BurstFrameIndices);
        
            % Save outputs
            mat_output_path = fullfile(root_folder, [output_file_name '.mat']);
            csv_output_path = fullfile(root_folder, [output_file_name '.csv']);
        
            save(mat_output_path, 'batch_network_bursts');
            writetable(batch_network_bursts, csv_output_path);
        end
    end
end