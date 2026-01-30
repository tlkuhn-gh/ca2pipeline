function network_burst_gui()
    addpath('NetworkBursts.m');
% NETWORK_BURST_GUI
% GUI for batch processing and visualization of network bursts

    %% Main figure
    fig = uifigure('Name', 'Network Burst Viewer', ...
                   'Position', [100 100 1100 600]);

    gl = uigridlayout(fig, [1 2]);
    gl.ColumnWidth = {300, '1x'};

        %% Left panel: controls
    ctrl = uigridlayout(gl, [13 1]);
    ctrl.RowHeight = { ...
        18, 22, ...   % Root folder
        18, 22, ...   % File prefix
        18, 22, ...   % Output file name
        18, 22, ...   % Min # neurons per burst
        18, 22, ...   % Burst frame range
        22, ...       % Include subfolders
        30, ...       % Run button
        '1x', ...     % File list
        22};          % Status label

    ctrl.Padding = [10 10 10 10];

    % Root folder
    uilabel(ctrl, 'Text', 'Path to FluoroSNNAP output folder');
    rootEdit = uieditfield(ctrl, 'text');

    % File prefix
    uilabel(ctrl, 'Text', 'Input file prefix (of all files w/calcium event data)');
    prefixEdit = uieditfield(ctrl, 'text', 'Value', 'analysis-');

    % Output file name (will save as .csv and .mat)
    uilabel(ctrl, 'Text', 'Output file name (will save as .csv and .mat)');
    outputFileEdit = uieditfield(ctrl, 'text', ...
        'Placeholder', 'Output file name', 'Value', 'network_bursts');

    % Percent neurons
    uilabel(ctrl, 'Text', 'Minimum # neurons for burst');
    minNeuronsEdit = uieditfield(ctrl, 'numeric', ...
        'Limits', [0 3000], ...
        'Value', 100);

    % Burst frame range
    uilabel(ctrl, 'Text', 'Burst frame range');
    frameRangeEdit = uieditfield(ctrl, 'numeric', ...
        'Limits', [0 Inf], ...
        'RoundFractionalValues', true, ...
        'Value', 10);

    % Include subfolders
    includeSubChk = uicheckbox(ctrl, ...
        'Text', 'Include subfolders');

    % Run button
    runBtn = uibutton(ctrl, 'push', ...
        'Text', 'Run batch processing', ...
        'ButtonPushedFcn', @runBatch);

    % File list
    fileList = uilistbox(ctrl, ...
        'ValueChangedFcn', @fileSelected);

    % Status label
    statusLabel = uilabel(ctrl, ...
        'Text', 'Idle');

    %% Right panel: plot
    ax = uiaxes(gl);
    title(ax, 'No file selected');
    xlabel(ax, 'Frame');
    ylabel(ax, 'Neuron Index');

    %% Shared state
    state.batch_results = [];
    fig.UserData = state;

    %% Callback: run batch processing
    function runBatch(~, ~)
        statusLabel.Text = 'Processing...';
        drawnow;

        try
            results = NetworkBursts.batch_process_network_bursts( ...
                rootEdit.Value, ...
                includeSubChk.Value, ...
                prefixEdit.Value, ...
                minNeuronsEdit.Value, ...
                frameRangeEdit.Value, ...
                outputFileEdit.Value);

            fig.UserData.batch_results = results;
            fileList.Items = results.FileName;

            statusLabel.Text = 'Done';
        catch ME
            statusLabel.Text = 'Error';
            uialert(fig, ME.message, 'Batch error');
        end
    end

    %% Callback: file selection
    function fileSelected(src, ~)
        idx = find(strcmp(src.Value, fig.UserData.batch_results.FileName));
        if isempty(idx)
            return;
        end

        entry = fig.UserData.batch_results(idx, :);

        S = load(entry.FullFilePath{1});
        data = S.data;

        cla(ax);
        hold(ax, 'on');

        % --- plot bursts first ---
        for i = 1:numel(entry.BurstFrameIndices{1})
            xline(ax, entry.BurstFrameIndices{1}(i), ...
                'r-', 'LineWidth', 2);
        end

        % --- plot raster ---
        num_neurons = numel(data.Spikes_cell);
        max_frame = 0;

        for n = 1:num_neurons
            frames = data.Spikes_cell{n};
            if ~isempty(frames)
                max_frame = max(max_frame, max(frames));
                plot(ax, frames, ...
                    repmat(n, size(frames)), ...
                    '.', 'Color', 'k', 'MarkerSize', 10);
            end
        end

        ylim(ax, [0.5 num_neurons + 0.5]);
        xlim(ax, [0 max_frame]);
        set(ax, 'YDir', 'reverse');

        title(ax, entry.FileName{1});
        hold(ax, 'off');
    end
end