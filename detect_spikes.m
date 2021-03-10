function [all_threshold_spiketimes,all_threshold_amplitudes] = detect_spikes(list_of_recordings, chunk_size, max_n_chunks, list_of_thresholds, resultpath, list_of_rec_resultpaths)



% total_n_channels = 0;
% n_channels_per_rec = [];
% n_channels_cumulative = [];

%prepare vector of channels that have invalid signals
%     invalid_channels_all_recordings = [];
%     path_snapshots_invalid = [resultpath '/invalid_channels/']
%     mkdir(path_snapshots_invalid)
%     path_low_MAD = [resultpath '/low_MAD_channels/']
%     mkdir(path_low_MAD)


%     for rec_idx = 1:length(list_of_recordings)
%         %recording_file = [datapath chip date recording_type recordings{rec_idx} filextension];
%         recording = mxw.fileManager(list_of_recordings{rec_idx});
%         [testdata, ~, ~] = recording.extractRawData(1, 10);
%         total_n_channels = total_n_channels + size(testdata,2);
%         n_channels_per_rec = [n_channels_per_rec size(testdata,2)];
%         n_channels_cumulative = [n_channels_cumulative total_n_channels];
%
%
%
%         %detect bad channels:
% %         %if the raw signal reaches <1 or >6340 over time,
% %         %it seems that there is something wrong with the
% %         %electrode
%         %nsamples = floor(recording.fileObj.dataLenSamples/20000)*20000
%         %n_test_snippets = 10;
%         snippet_size = 10000;
%         snipped_idxs = 5;
%         snippet_start_points = (101 * 20000) + snipped_idxs.*snippet_size   %randi(nsamples-snippet_size,n_test_snippets,1);
%         path_snapshots = '/home/euler/MAXONE/results_May18/May8/test/'
%         mkdir(path_snapshots);
%         for j = 1:length(snippet_start_points)%size(snippet_start_points,1)
%             rawDataPart = recording.extractRawData(snippet_start_points(j)+1000+120, snippet_size/50/5);
% %            for e = 1:size(rawDataPart,2)
%                 %if (find(rawDataPart(:,e) < 1 |  rawDataPart(:,e) >= 6340))
% %                     invalid_channels_all_recordings = [invalid_channels_all_recordings e];
% %                     invalid_channels_this_recording = [invalid_channels_this_recording e];
%             f = figure();
%             plot(rawDataPart);
%             pictureName = ['raw_data_snippet_' int2str(snipped_idxs(j))];
%             savepng( 'Directory', path_snapshots, 'FileName' , pictureName );
%             hold off;
%             close(f);
%                 %end
% %            end
%         end
%        dlmwrite([list_of_rec_resultpaths{rec_idx} '/invalid_indices.csv'],unique(invalid_channels_this_recording));
%     end
%

%%create and prepare data structures

    %determine the total number of recorded channels
    recording_file = list_of_recordings{1}
    recording = mxw.fileManager(recording_file);
    [data, ~, ~] = recording.extractBPFData(1, 1000); %extract chunks of bandpass filtered data
    total_n_channels = size(data, 2)
    
    all_threshold_spiketimes = {};
    all_threshold_amplitudes = {};

%     total_electrodes = zeros(total_n_channels, 1);
%     total_x = zeros(total_n_channels, 1);
%     total_y = zeros(total_n_channels, 1);
%     total_average_amplitudes = zeros(total_n_channels, 1);
%     total_rec_indices =  zeros(total_n_channels, 1);


    for threshold = list_of_thresholds
        
        %create and prepare the sparse structure to save the spiketimes and
        %corresponding amplitudes
        total_spiketimes = {};
        total_amplitudes = {};
        for i = 1:total_n_channels
            total_spiketimes{end+1} = []; %these empty entries will be filled with lists of spiketimes
            total_amplitudes{end+1} = []; %these empty entries will be filled with lists of corresponding amplitudes
        end
        
        
        all_threshold_spiketimes{end+1} = total_spiketimes;
        all_threshold_amplitudes{end+1} = total_amplitudes;
        
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% load data in chunks, filter, determine spike times, and save in sparse structure %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load the (filtered) data chunk per chunk; determine the spike times and fill the sparse structure

for rec_idx = 1:length(list_of_recordings)
    
    %recording_file = [datapath chip date recording_type recordings{rec_idx} filextension];
    recording_file = list_of_recordings{rec_idx}
    recording = mxw.fileManager(recording_file);
    map = recording.rawMap.map;
    
    % to collect the invalid channels indices
    invalid_channels_this_recording = [];
    
    %         %save the meta-information
    %         if rec_idx == 1
    %             range_in_total = 1 : n_channels_cumulative(1);
    %         else
    %             range_in_total = n_channels_cumulative(rec_idx-1) + 1 : n_channels_cumulative(rec_idx);
    %         end
    %         total_electrodes(range_in_total, 1)  = map.electrode;
    %         total_x(range_in_total, 1) = map.x;
    %         total_y(range_in_total, 1) = map.y;
    %         total_rec_indices(range_in_total, 1) = repelem(rec_idx, length(range_in_total)) ;
    
    effective_n_chunks = 0;
    %detect the spikes by loading the data chunk by chunk
    for c = 1:max_n_chunks
        try
            %extract filtered data
            [data, ~, ~] = recording.extractBPFData(1 + (c-1)*chunk_size, chunk_size); %extract chunks of bandpass filtered data
            [rawData, ~, ~] = recording.extractRawData(1 + (c-1)*chunk_size, chunk_size); %extract chuncks of raw data to detect bad sensors
            
            %compute MAD for all electrodes
            median_abs_dev = mad(data, 1);
            
            % determine spike times
            for e = 1:size(data,2)%n_channels_per_rec(rec_idx) %for all electrodes
                %electrode_index_in_total =  range_in_total(e);
                
                if isempty(invalid_channels_this_recording) | isempty(find(invalid_channels_this_recording == e)) %the channel was not found to be invalid based on the previous chunks
                    
                    for threshold_idx = 1:length(list_of_thresholds) %for all thresholds
                        
                        threshold =  list_of_thresholds(threshold_idx);
                        
                        for i = 3:chunk_size-3 %for all timesteps in this chunk
                            
                            if (rawData(i,e)<1 | rawData(i,e)>= 6340) %the channel is invalid: add to the list of invalid channels and stop the processing
                                invalid_channels_this_recording = [invalid_channels_this_recording e];
                                break;
                            elseif (data(i,e) < -threshold * median_abs_dev(e) & data(i,e) < data(i-1,e) & data(i,e) < data(i-2,e) & data(i,e) < data(i+1,e) & data(i,e) < data(i+2,e) )
                                %all_threshold_spiketimes{threshold_idx}{electrode_index_in_total}(end+1) = i + (c-1)*chunk_size;
                                %all_threshold_amplitudes{threshold_idx}{electrode_index_in_total}(end+1) = data(i,e);
                                all_threshold_spiketimes{threshold_idx}{e}(end+1) = i + (c-1)*chunk_size;
                                all_threshold_amplitudes{threshold_idx}{e}(end+1) = data(i,e);
                                
                            end
                        end
                    end
                end
                
                if ~isempty(invalid_channels_this_recording) &  find(invalid_channels_this_recording == e) %the channel is invalid
                    for threshold_idx = 1:length(list_of_thresholds)
                        %all_threshold_spiketimes{threshold_idx}{electrode_index_in_total} = [];
                        %all_threshold_amplitudes{threshold_idx}{electrode_index_in_total} = [];
                        all_threshold_spiketimes{threshold_idx}{e} = [];
                        all_threshold_amplitudes{threshold_idx}{e} = [];
                    end
                end
            end
            effective_n_chunks = effective_n_chunks + 1
            clearvars data
            clearvars rawData
        end
    end
    dlmwrite([list_of_rec_resultpaths{rec_idx} '/invalid_indices.csv'],unique(invalid_channels_this_recording));
    
    %%%%%%% remove glitches
    [~, largest_threshold_idx] = max(list_of_thresholds);
    if ~isempty(all_threshold_spiketimes{largest_threshold_idx})
        
        disp('detecting and removing recording glitches...')
        
        %time_intervals_glitches = detect_glitches(all_threshold_spiketimes{largest_threshold_idx}(range_in_total), all_threshold_amplitudes{largest_threshold_idx}(range_in_total), map);
        
        time_intervals_glitches = detect_glitches(all_threshold_spiketimes{largest_threshold_idx}, all_threshold_amplitudes{largest_threshold_idx}, map);
        dlmwrite([list_of_rec_resultpaths{rec_idx} '/time_intervals_glitches.csv'],cell2mat(time_intervals_glitches));
        
        
        if ~isempty(time_intervals_glitches)
            if ~isempty(time_intervals_glitches{1})%to do: fix the output of time_intervals_glitches to be uniform for all cases
                for glitch_interval_idx = 1:length(time_intervals_glitches)
                    glitch_interval = time_intervals_glitches{glitch_interval_idx}
                    %spiketimes = all_threshold_spiketimes{largest_threshold_idx}(range_in_total);
                    spiketimes_current_threshold = all_threshold_spiketimes{largest_threshold_idx};
                    
                    %remove all the spiketimes that happen within the glitch interval
                    %for all electrodes in the recording
                    for e = 1:size(spiketimes_current_threshold, 2)
                        spiketimes_electrode = spiketimes_current_threshold{e}./20000;
                        correct_index = e;%range_in_total(e);
                        spiketimes_to_remove = [];
                        for s = 1:length(spiketimes_electrode)
                            if (spiketimes_electrode(s) >= glitch_interval(1) & spiketimes_electrode(s) <= glitch_interval(2))
                                spiketimes_to_remove = [spiketimes_to_remove s];
                            end
                        end
                        all_threshold_spiketimes{largest_threshold_idx}{correct_index}(spiketimes_to_remove) = [];
                        all_threshold_amplitudes{largest_threshold_idx}{correct_index}(spiketimes_to_remove) = [];
                    end
                end
            end
        end
        
    else
        disp('no spikes detected at largest threshold, no need to remove glitches')
    end
    
    
    number_of_chunks_in_recording = effective_n_chunks
    number_of_seconds_in_recording = effective_n_chunks * chunk_size / 20000
    
    rec_csv = [resultpath 'rec.csv'];
    if exist(rec_csv, 'file') == 2
        system(['rm -rf ' rec_csv])
    end
    dlmwrite(rec_csv,[number_of_chunks_in_recording, number_of_seconds_in_recording]);
    

    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% write the results to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%% write the results for each threshold to file
    
    %
    %     for threshold_idx = 1:length(list_of_thresholds) %for all thresholds
    %
    %         total_spiketimes = all_threshold_spiketimes{threshold_idx};
    %         total_amplitudes = all_threshold_amplitudes{threshold_idx};
    %
    %         threshold = list_of_thresholds(threshold_idx);
    %
    %         spiketimesfile_csv = [resultpath 'sensor_spiketimes_'  int2str(threshold)  '.csv']
    %         amplitudesfile_csv =  [resultpath 'sensor_amplitudes_'  int2str(threshold)  '.csv']
    %         %since we append to this file, check whether it should be
    %         %removed -> so we don't append to an older version of the file
    %         if exist(spiketimesfile_csv, 'file') == 2
    %             system(['rm -rf ' spiketimesfile_csv])
    %         end
    %
    %         if exist(amplitudesfile_csv, 'file') == 2
    %             system(['rm -rf ' amplitudesfile_csv])
    %         end
    %
    %         nspikes = zeros(length(total_spiketimes), 1);
    %         mean_ampl =  zeros(length(total_spiketimes), 1);
    %         nspikes_csv = [resultpath 'sensor_nspikes_'  int2str(threshold)  '.csv'];
    %         meanAmpl_csv = [resultpath 'sensor_meanAmpl_'  int2str(threshold)  '.csv'];
    %         %since we append to this file, check whether it should be
    %         %removed -> so we don't append to an older version of the file
    %         if exist(nspikes_csv, 'file') == 2
    %             system(['rm -rf ' nspikes_csv]);
    %         end
    %
    %         if exist(meanAmpl_csv, 'file') == 2
    %             system(['rm -rf ' meanAmpl_csv]);
    %         end
    %
    %
    %         for k = 1:length(total_spiketimes)
    %            nspikes(k,1) = length(total_spiketimes{k});
    %            mean_ampl(k,1)  = mean(total_amplitudes{k});
    %            dlmwrite(nspikes_csv, length(total_spiketimes{k}), '-append','delimiter',',');
    %
    %            if(~isempty(total_spiketimes{k}))
    %               dlmwrite(meanAmpl_csv, mean(total_amplitudes{k}), '-append','delimiter',',');
    %            else
    %               dlmwrite(meanAmpl_csv, 0, '-append','delimiter',',');
    %            end
    %
    %            if isempty(total_spiketimes{k})
    %               dlmwrite(spiketimesfile_csv, ' ', '-append','delimiter',',');
    %               dlmwrite(amplitudesfile_csv, ' ', '-append','delimiter',',');
    %            else
    %               dlmwrite(amplitudesfile_csv,total_amplitudes{k},'-append','delimiter',',', 'precision', 7);
    %               dlmwrite(spiketimesfile_csv,total_spiketimes{k},'-append','delimiter',',', 'precision', 7);
    %            end
    %         end
    %
    %         figure()
    %         histogram(nspikes/number_of_seconds_in_recording);
    %         title(['Detected Spike Frequency, with threshold ' int2str(threshold)])
    %         xlabel('Frequency (Hz)')
    %         ylabel('Number of Sensors')
    %         pictureName = ['sensor_spikes_histogram_' int2str(threshold)];
    %         savepng( 'Directory', resultpath, 'FileName' , pictureName );
    %         hold off
    %
    %         figure()
    %         histogram(abs(mean_ampl));
    %         title(['Detected mean spike amplitude, with threshold ' int2str(threshold)])
    %         xlabel('Abs. Value Amplitude (microV)')
    %         ylabel('Number of Sensors')
    %         pictureName = ['sensor_amplitudes_histogram_' int2str(threshold)];
    %         savepng( 'Directory', resultpath, 'FileName' , pictureName );
    %         hold off
    %     end
    %
    %

    
    
    
    
    
    
    
    
    
end

