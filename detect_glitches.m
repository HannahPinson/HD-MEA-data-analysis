function time_intervals_glitches = detect_glitches(spiketimes, amplitudes, map)
   
    % a glitch is here detected as a sudden, unphysical increase in
    % activity for a very short period of time that happens at the same
    % over the whole recording (thus in both squares, in case of a two
    % square recording)
    
    
    %divide the recording in two squares (can be done even if it is a
    %single square recording)
    rel_middle = (max(map.x) - min(map.y))/2;
    first_square = map.x< (min(map.x) + rel_middle); %indices in map correspoding to first square
    second_square = map.x >= (min(map.x) + rel_middle);%indices in map correspoding to second square
   
    spiketimes_sets = {}; 
    spiketimes_sets{end+1} =spiketimes(first_square);
    spiketimes_sets{end+1} =spiketimes(second_square);
    
    amplitudes_sets = {}; 
    amplitudes_sets{end+1} = amplitudes(first_square);
    amplitudes_sets{end+1} = amplitudes(second_square); 

    binning_factor = 50 * 20; %number of samples
    sampling_frequency = 20000;
    
    extreme_value_bins_per_square = {};
    
    for i = 1:2
        spiketimes_set = spiketimes_sets{i};
        amplitudes_set = amplitudes_sets{i};
        merged_spiketimes_list = []; 
        merged_amplitudes_list = [];
        for c = 1:length(spiketimes_set)
            channel_spiketimes = spiketimes_set{c};
            channel_amplitudes = amplitudes_set{c};
            merged_spiketimes_list = [merged_spiketimes_list channel_spiketimes];
            merged_amplitudes_list = [merged_amplitudes_list channel_amplitudes];
        end
        

        histo_spiketimes = histogram(merged_spiketimes_list, 'BinWidth',binning_factor);
        activity_values = histo_spiketimes.Values;
        activity_bins = histo_spiketimes.BinEdges;
        
        [sorted_spiketimes, sorting_indices] = sort(merged_spiketimes_list);
        sorted_amplitudes = merged_amplitudes_list(sorting_indices); 
        
        median_abs_dev_activity = mad(activity_values);
        median_abs_dev_amplitudes = mad(merged_amplitudes_list);
        extreme_value_bins = [];
        
        for index = 1:length(activity_values)
            act_value = activity_values(index);
            bin = activity_bins(index):activity_bins(index+1);
            corresponding_idx = ismember(sorted_spiketimes, bin);
            ampl_value = max(abs(sorted_amplitudes(corresponding_idx)));
%             if size(find(corresponding_idx==1),2)~=0
%                 find(corresponding_idx==1)
%                 disp("succes")
%                 ampl_value
%             end
%             

            if (act_value > 3 * median_abs_dev_activity & abs(ampl_value) > abs( 10 * median_abs_dev_amplitudes) )
                extreme_value_bins = [extreme_value_bins activity_bins(index)/sampling_frequency];
            end
        end

        extreme_value_bins_per_square{end+1} = extreme_value_bins;
    end
    
    
    
    shared_extreme_bins = intersect(extreme_value_bins_per_square{1}, extreme_value_bins_per_square{2});
    time_intervals_glitches = {};
    for glitch = shared_extreme_bins
        time_intervals_glitches{end+1} = [glitch glitch+binning_factor/sampling_frequency];
    end
        
end

