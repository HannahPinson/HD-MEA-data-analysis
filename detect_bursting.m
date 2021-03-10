%% detect bursting 

function [all_within_burst_intervals, all_between_burst_intervals, all_burst_durations, all_n_spikes_per_burst, all_n_bursts, burst_duration_per_sensor, n_bursts_per_sensor, mean_n_spikes_per_sensor] = detect_bursting(spiketimes)

    sampling_frequency = 20000;

    all_within_burst_intervals = [];
    all_between_burst_intervals = []; 
    all_burst_durations = [];
    all_n_spikes_per_burst = []; 
    all_n_bursts = [];
    
    n_bursts_per_sensor = []; 
    burst_duration_per_sensor = []; 
    mean_n_spikes_per_sensor = []; 
    
    for s = 1:length(spiketimes)
        
        selected_spiketimes = spiketimes{s};%spiketimes{392};
       

        if length(selected_spiketimes) > 5

            ISI = diff(selected_spiketimes)./sampling_frequency;
            histBinWidth = 0.005;

            %figure()
            %histogram(ISI)
    %         
    %         histo = histogram(ISI, 'BinWidth',histBinWidth);
    %         binned_ISI_counts = histo.Values; 
    %         binned_ISI_edges = histo.BinEdges;
            %hold off

            [binned_ISI_counts,binned_ISI_edges] = histcounts(ISI,'BinWidth',histBinWidth);




            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %%%%%% determine detection threshold %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %calculate cumulative moving average 

            CMA = zeros(size(binned_ISI_counts));
            sum = 0;
            for i = 1:length(CMA)
                sum = sum + binned_ISI_counts(i); 
                CMA(i) = sum / i; 
            end

            [max_CMA, idx_max_CMA] = max(CMA);
            x_m = binned_ISI_edges(idx_max_CMA) + histBinWidth/2;

            %calculate ISI threshold
            %the ISI threshold xt, xt > xm, for burst detection is found at the mid time point of the ISI bin for which the value of the CMA curve is the closest to ?·CMAm

            skew = skewness(ISI);
            
            alpha = 1;
            if skew > 1 & skew < 4
                alpha = 0.7;
            elseif  skew > 4 & skew < 9
                alpha = 0.5;
            elseif skew > 9
                alpha = 0.3;
            end
            
            alpha2 = 0.5;
            if skew > 4 & skew < 9
                alpha2 = 0.3;
            elseif skew > 9
                alpha2 = 0.1;
            end

            y_t = alpha * max_CMA; %CMA value equal to alpha times the max of the CMA function
            %now we have to look for the corresponding bin (=x_t) on the right of the
            %bin of the CMA max (=x_m). We assume CMA is, after the max value is
            %reached, monotonically decreasing. 
            bin_index = 1;
            for i = idx_max_CMA+1:length(CMA)-1
                if(CMA(i-1) >= y_t) & (CMA(i+1) <= y_t)
                    bin_index = i;
                    break;
                end

            end
            
            index_x_t = bin_index; 
            x_t = binned_ISI_edges(bin_index) + histBinWidth/2;
            
            
            %same procedure for the second threshold
            y_t2 = alpha2 * max_CMA; %CMA value equal to alpha times the max of the CMA function
            for i = idx_max_CMA+1:length(CMA)-1
                if(CMA(i-1) >= y_t2) & (CMA(i+1) <= y_t2)
                    bin_index = i;
                    break;
                end
            end

            index_x_t2 = bin_index; 
            x_t2 = binned_ISI_edges(bin_index) + histBinWidth/2;
            



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %%%%%%%%%%%%%%%% burst detection %%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            burst_spikes = {};
            non_burst_spikes = []; %keep track of the spiketimes that do not belong to bursts for use in correction step
            burst_intervals = {}; 

            spikes_in_burst = 1; 
            current_burst_spikes = [selected_spiketimes(1)/sampling_frequency];

            %detect the bursts
            for i = 1:length(selected_spiketimes)-1
                current_spiketime = selected_spiketimes(i)/sampling_frequency; 
                next_spiketime = selected_spiketimes(i+1)/sampling_frequency;

                %inter spike interval smaller than the threshold: within burst
                if (next_spiketime - current_spiketime <= x_t)
                    if (spikes_in_burst == 1) %current spike is the first spike of potential burst
                        potential_start_burst = current_spiketime; 
                    end
                    spikes_in_burst = spikes_in_burst+1;
                    current_burst_spikes = [current_burst_spikes next_spiketime]; %the next spike is part of this burst

                %the next spike is further away than the threshold: end of burst, or no burst at all
                else 
                    if (spikes_in_burst >= 3) %it was a burst, but the current spike is the last one
                        burst_intervals{end+1} = [potential_start_burst current_spiketime];
                        burst_spikes{end+1} = current_burst_spikes;
                    else
                        non_burst_spikes = [non_burst_spikes current_burst_spikes];
                    end
                    %the next spike could again be the start of a burst
                    spikes_in_burst = 1; %start counting again from 1
                    current_burst_spikes = [next_spiketime];
                      
                end 
            end
            
            
            %%correction step
            
            %add all burst related spikes (non burst spikes within second threshold
            %of a burst) to corresponding bursts
            
            
            burst_related_spikes = [];
            
%             for i = 1:length(non_burst_spikes)-1
%                 
%                 non_burst_spike = non_burst_spikes(i);
%                 
%                 for j = 1:length(burst_intervals)
%                     
%                     burst_interval = burst_intervals{j};
%                     if burst_interval(1) >= non_burst_spike & abs(burst_interval(1) -  non_burst_spike) <= x_t2 %the burst related spike should be added in front of the burst             
%                         burst_spikes{j} = [non_burst_spike burst_spikes{j}];
%                         burst_intervals{j}(1) = non_burst_spike; 
%                         burst_related_spikes = [burst_related_spikes non_burst_spike];
%                     elseif non_burst_spike >= burst_interval(2) & abs(non_burst_spike - burst_interval(2)) <= x_t2 %the burst related spike should be added at the end of the burst
%                         burst_spikes{j} = [burst_spikes{j} non_burst_spike ];
%                         burst_intervals{j}(2) = non_burst_spike; 
%                         burst_related_spikes = [burst_related_spikes non_burst_spike];
%                     end
%                 end
%             end
%             
%             
%             %merge all bursts within second threshold from one another 
%             
            indices_of_bursts_to_remove = [];
            
            for j = 1:length(burst_intervals)-1

                current_burst_end = burst_intervals{j}(2);
                next_burst_begin =  burst_intervals{j+1}(1);
                next_burst_end = burst_intervals{j+1}(2);
                if abs(next_burst_begin - current_burst_end) <= x_t2
                    %merge bursts
                    burst_intervals{j}(2) =  next_burst_end;
                    burst_spikes{j}  = [burst_spikes{j}  burst_spikes{j+1}];
                    
                    %keep track of the burst that needs to be removed
                    indices_of_bursts_to_remove = [indices_of_bursts_to_remove j+1];
                    
                end
            end
            
            
            shift = 0; 
            for k = indices_of_bursts_to_remove
                burst_intervals(k-shift) = [];
                burst_spikes(k-shift) = [];
                shift = shift+1;
            end
            
            
            
            

            % plot a single spike train
            % 
% 
%             if (mod(s, 10) == 0)%length(indices_of_bursts_to_remove)>1%
%                 figure()
%                 endtime = 60 * 20000;
%                 time_range= 1:endtime;
%                 plot(ismember(time_range, selected_spiketimes), 'k')
%                 hold on 
%                                 
%                 plot(ismember(time_range, burst_related_spikes.*sampling_frequency), 'r');
%                 hold on
%     
%                 burst_timepoints = [];
%                 for i = 1:length(burst_intervals)
%                     burst_interval_indices = burst_intervals{i}.*sampling_frequency;
%                     burst_range_indices = burst_interval_indices(1):burst_interval_indices(2);
%                     burst_timepoints = [burst_timepoints burst_range_indices];
%                     %burst_line(burst_interval_indices(1):burst_interval_indices(2)) = ones(size(burst_range_indices)).*1.2; 
%                 end
%     
%                 scatter(burst_timepoints, ones(size(burst_timepoints)).*1.2, 10, 'filled', 'b'); 
%                 ylim([0 1.5])
%                 xlim([1 endtime]);
%                 hold off
%                 
%                 pictureName = ['merged_burst_detection_', int2str(s)];
%                 savepng( 'Directory', '/Volumes/Hannah_Maxwell_Data/Hannah_Maria_data_010719/resultsMaria/testBurstDetection', 'FileName' , pictureName );
%                 
%              end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %%%%%%%%%%%%% compute burst metrics %%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %compute all metrics

            within_burst_intervals = [];
            between_burst_intervals = []; 
            burst_durations = [];
            n_spikes_per_burst = []; 
            n_bursts = length(burst_spikes);

            for i = 2:n_bursts

                within_burst_intervals = [within_burst_intervals diff(burst_spikes{i})];
                if i<length(burst_spikes)-1
                    between_burst_intervals = [between_burst_intervals (burst_intervals{i+1}(1) - burst_intervals{i}(1))];
                end

                burst_durations = [burst_durations (burst_intervals{i}(2) -  burst_intervals{i}(1))];
                n_spikes_per_burst = [n_spikes_per_burst length(burst_spikes{i})];

            end

            all_within_burst_intervals = [all_within_burst_intervals within_burst_intervals];   
            all_between_burst_intervals = [all_between_burst_intervals between_burst_intervals]; 
            all_burst_durations = [ all_burst_durations burst_durations];
            all_n_spikes_per_burst = [ all_n_spikes_per_burst n_spikes_per_burst]; 
            all_n_bursts = [all_n_bursts n_bursts];
            n_bursts_per_sensor = [n_bursts_per_sensor n_bursts];
            
            if n_bursts < 2
                burst_duration_per_sensor = [ burst_duration_per_sensor 0];
                mean_n_spikes_per_sensor = [mean_n_spikes_per_sensor 0]; 
            else
                burst_duration_per_sensor = [ burst_duration_per_sensor mean(burst_durations)];
                mean_n_spikes_per_sensor = [mean_n_spikes_per_sensor mean(n_spikes_per_burst)]; 
            end

        else
            n_bursts_per_sensor = [n_bursts_per_sensor 0]; 
            burst_duration_per_sensor = [ burst_duration_per_sensor 0]; 
            mean_n_spikes_per_sensor = [mean_n_spikes_per_sensor 0];
        end
    end
    
  


end













%




% 
% C_matrix_path = [folder 'CMatrices_binned_10ms_1000ms/'];
% order = 100;
% C_total = readCMatrices(C_matrix_path, order);
% nvar = size(C_total, 1)
% somaNumbers = 1:nvar;
% 
% significant_pairs = {};
% i = selected_index %1:length(somaNumbers)
% somaNumber = somaNumbers(i);
% j = i
% otherSomaNumber = somaNumbers(j);
% normalization = sqrt(C_total(i,i) * C_total(j,j));
% 
% %read out the cross correlation values 
% auto_corr =  zeros([1,order]);
% for o = 1:order
%     auto_corr(1,o) = C_total(somaNumber, nvar * (o-1) + otherSomaNumber)/normalization;
% end
% 
% figure()
% plot(auto_corr(1, 2:end));


% 
% options = statset('MaxIter',1000);
% f = fitgmdist(transpose(auto_corr(1, 2:end)),1,'RegularizationValue',0.01, 'Options', options);
% 
% 
% %xgrid = linspace(2,1,100)';
% xgrid = -0.1:0.001:0.1;
% xgrid = transpose(xgrid);
% figure() 
% n1 = makedist('normal',f.mu(1),sqrt(f.Sigma(1)));
% p = f.ComponentProportion;
% y_1 = p(1)*pdf(n1,xgrid); 
% plot(xgrid,y_1,'r-','LineWidth', 2); 
%         



