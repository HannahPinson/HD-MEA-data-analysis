addpath(genpath(pwd))
warning('off','all')

filextension = '.raw.h5';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% parameters of the dataset %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%March 18 data
% chips = {'2057_150k', '2058_200k', '2060_75k', '2061_50k'}
% dates = {'180326','180328','180330', '180402','180404'} %{'180316', '180319', '180321', '180326'}
% recNumberFormat = '%03d';
% basename = ''
% squares = [1]
% type = ''
% recordings = 0:5;

% % %May 18 data
% chips = {'2170_LAM77-1'}%{'2057_LAM78', '2169_LAM78', '2061_LAM77', '2170_LAM77-1', '2199_LAM77-1', '2157_LAM77'}
% dates = {'062018'} %{'051618', '052118', '052318', '052518', '052818'}
% recNumberFormat = '%01d'%
% recordings = 0:24%[0:13 15:16 23:24] %%0:19%0:24%strsplit(int2str(0:5)
% basename = 'denseConfig_'%'2199_LAM77-1_network' %''
% squares = [1,2] %[1]%
% type = '23X22'; %'network_traces_10min_300Hz'%



% chips = {'2170_LAM77-1'}%{'2199_LAM77-1'}% {'2169_LAM78'}%  '2061_LAM77', '2199_LAM77-1'}
% dates = {'062018'} %{'051618', '052118', '052318', '052518', '052818'}
% recNumberFormat = '%01d' %network '%01d'% dense: '%03d';
% recordings = 0%0:19%[0:13 15:16 23:24] %%0:19%0:24%strsplit(int2str(0:5)
% basename = 'network_session0'%'2199_LAM77-1_network' %'2169_LAM78_network'%'2199_LAM77-1_network' %''
% squares = [1] %[1,2] %%
% type = 'network' %'network_traces_10min_300Hz' %'20x20_2min'; %%


%%january 19 data
% chips = {'A'}
% dates = {'190121'}%{'190121'}%
% recNumberFormat = '%01d'
% recordings = 9 %0:24
% basename = ''
% squares = [1, 2]


%%may 19 full dense data
% recNumberFormat = '%01d'
% recordings = 0:32
% basename = ''
% squares = [1 2]
% chips = {'100'}
% dates = {'190619' }%{'190602', '190605', '190607', '190610', '190611', '190612' }
% type = 'spontaneous_before_training'


%%may 19 parameter testing data
% recNumberFormat = '%01d'
% recordings = 0
% basename = ''
% squares = [1]
% chips = {'50'}
% dates = {'190530'}
% type = 'paremeter_testing'


%%data Maria

% chips = {'2454', '2467', '2468'};
% group = 'group1_2454_2467_2468';
% dates = {'080618','080618','080918','081318','081618','082018','082318', '082718', '083018'}
%
% chips = {'2455', '2465', '2469'}
% group = 'group2_2455_2465_2469'
% dates = {'081518'}; %{'082118', '082418','082818','083118','090418'}; %{'081518', '082118', '082418','082818','083118','090418'};
% %
% recNumberFormat = '%01d';
% recordings = 0:24;
% basename = 'denseConfig_';
% type='23X22';
% squares = [1,2];


%%data Maria part 2

% group = 'group3_3722_3724_3726'
% chips = {'3722', '3724', '3726'}
% dates = {'190404', '190408', '190409', '190415'}%{'190408', '190409', '190415', '190418'} %{'190321', '190325', '190328', '190401', '190415'}
% types_to_try =  {'Activity_Scan'}%{'23X22','23x22', '20x20', '20X20', 'Activity_Scan'}
% ids_to_try = {'000068', '000063', '000059', '000126', '000113', '000133', '000228', '000216', '000222'}%{'000068', '000063', '000059','000117','000135', '000129','000219','000266', '000259'}

%  group = 'group5_3723_3725'
%  chips = {'3723', '3725'}
%  dates = {'190416'}
%  types_to_try = {'20x20', '23X22', '23x22', '23x23', '23X23', 'Activity_Scan'}
%  ids_to_try = {'000249', '000241'}%{'000198', '000192'}%{'000045'}
%  type = 'dense'

%group = 'group4_2469'
% chips = {'2469'}
% dates = {'190410'}%{'190324'}%{'190413'}
% type = 'dense'
% types_to_try = {'20x20', '23X22', '23x22', '23x23', '23X23', 'Activity_Scan'}
% ids_to_try = {'000097'}%{'000045'}


%data Maria network
%
% group = 'network'
% chips = {'2454', '2467', '2468', '2455', '2465', '2469'}%{'3723', '3725'} %{'2469'}
% dates = {''}%{'190412'} %{'190410'}
% type = 'network'
% recNumberFormat = -1%'%01d';
% recordings = 0;
% basename = -1;'recording';
% squares = [1];
%


% array_of_matching_paths = findfiles('Fullscan', '/home/euler/DATA1/DataMaria_part2/190401/')
% chips = {'3722', '3724','3726'}%{'3722','3724','3726'}
% dates = {'190401'}
%
% for l = array_of_matching_paths
%     l
% end
%
% global_path = ['/home/euler/MAXONE/resultsMaria/']
% experiment_group3 = 'group3_3722_3724_3726/resultsJune2/'



rec_type = 'axon_recording'%'23x22_scan'%light_sliding_each_10s'%'23x22_scan' %'mobile_phone_alternating_2'
array_of_matching_paths = findfiles(rec_type, '/home/euler/DATA2/levin_oct19/')
chips =  {'4163'} %{'4163'} %{'4162','4163','4164', '4165' }%;%{'4164','4162' }%
dates = {'191011'}%{'191012'}%,'191014'}%{'191008','191009', '191010'} % {'191011'}

for l = array_of_matching_paths
    l
end

global_path = ['/home/euler/MAXONE/resultsLevin/']
mkdir(global_path)


%
% rec_type = '190121A1255'%light_sliding_each_10s'%'23x22_scan' %'mobile_phone_alternating_2'
% array_of_matching_paths = findfiles(rec_type, '/home/euler/DATA2/MEA/')
% rec_type = '23x22_scan'
% chips =  {'A'} %{'4162','4163','4164', '4165' }%;%{'4164','4162' }%
% dates = {'190121'}%,'191014'}%{'191008','191009', '191010'} % {'191011'}
%
% for l = array_of_matching_paths
%     l
% end
%
% global_path = ['/home/euler/MAXONE/resultsJanuary/']
% mkdir(global_path)


%
%keep only the files
array_of_matching_files = array_of_matching_paths(contains(array_of_matching_paths,'raw.h5'))
recNumberFormat = '%03d'

if strcmp(rec_type,'23x22_scan') | strcmp(rec_type,'block_scan')
    squares = [1,2];
else
    squares = [1];
end

for l = array_of_matching_files
    l
end


%for debugging, do not use all files
% array_of_matching_files = array_of_matching_files(contains(array_of_matching_files, ['/' chips{1} '/']))
% array_of_matching_files = array_of_matching_files(contains(array_of_matching_files, ['/' dates{1} '/']))
% array_of_matching_files = array_of_matching_files(1:2)
% 
% for l = array_of_matching_files
%     l
% end






%data May 18 network
%
% group = 'network'
% chips = { '2057_LAM78', '2061_LAM77', '2157_LAM77',	'2169_LAM78', '2170_LAM77-1', '2199_LAM77-1'}%{'3723', '3725'} %{'2469'}
% dates = {''}%{'190412'} %{'190410'}
% type = 'network'
% recNumberFormat = -1%'%01d';
% recordings = 0;
% basename = -1;'recording';
% squares = [1];
%
% array_of_matching_paths = findfiles('network', '/home/euler/DATA4/MEA_rawDATA/June2018/EG180412/')
%
% %keep only the files
% array_of_matching_files = array_of_matching_paths(contains(array_of_matching_paths,'.raw.h5'))




% %recordings = -1 %defined later, depending on type of recording
% %basename = -1 %defined later, depending on type of recording %'denseConfig_';
% %recNumberFormat =  -1 %defined later, depending on type of recording '%01d';
% type = 'fullscan'
% %squares = -1 %defined later, depending on type of recording; [1,2];

% group = 'group4_2469'
% chips = {'2469'}
% dates = {'031219', '190318'}%{'190408', '190409', '190415', '190418'} %{'190321', '190325', '190328', '190401', '190415'}
% types_to_try =  {'23X22'}%{'23X22','23x22', '20x20', '20X20', 'Activity_Scan'}
% ids_to_try = {}%{'000068', '000063', '000059','000117','000135', '000129','000219','000266', '000259'}
% %recordings = -1 %defined later, depending on type of recording
% %basename = -1 %defined later, depending on type of recording %'denseConfig_';
% %recNumberFormat =  -1 %defined later, depending on type of recording '%01d';
% type = '23X22'
% %squares = -1 %defined later, depending on type of recording; [1,2];

% group = 'group5_2469'
% chips = {'2469'}
% dates = {'031219', '190318'}%{'190408', '190409', '190415', '190418'} %{'190321', '190325', '190328', '190401', '190415'}
% types_to_try =  {'23X22'}%{'23X22','23x22', '20x20', '20X20', 'Activity_Scan'}
% ids_to_try = {}%{'000068', '000063', '000059','000117','000135', '000129','000219','000266', '000259'}
% %recordings = -1 %defined later, depending on type of recording
% %basename = -1 %defined later, depending on type of recording %'denseConfig_';
% %recNumberFormat =  -1 %defined later, depending on type of recording '%01d';
% type = '23X22'
% %squares = -1 %defined later, depending on type of recording; [1,2];


%all result paths dopaminergic

%
% global_path = ['/home/euler/MAXONE/resultsMaria/']
%
%
% experiment_group1 = 'group1_2454_2467_2468/resultsApril29/'
% experiment_group2 = 'group2_2455_2465_2469/resultsApril29/'
% experiment_group3 = 'group3_3722_3724_3726/resultsJune2/'
% experiment_group4 = 'group4_2469/resultsJune2/'
% experiment_group5 = 'group5_3723_3725/resultsJune2/'
%
% basepaths = {experiment_group1,experiment_group1,experiment_group1, experiment_group2, experiment_group2, experiment_group2, experiment_group3,experiment_group3,experiment_group3, experiment_group4, experiment_group5,experiment_group5}



%
% data_to_read = {
%     {[experiment_group3  '190401/' '3722/' 'dense/'],[experiment_group3 '190401/' '3724/' 'dense/'],[experiment_group3  '190401/' '3726/' 'dense/']  }
%     }
% rec_types = {{'2squares','2squares','2squares'}}
%
% basepaths = {experiment_group3}

% data_to_read = {
%
% {[experiment_group1 '080618/' '2454/' '23X22/' ],
% [experiment_group1 '080918/' '2454/' '23X22/' ],
% [experiment_group1 '081318/' '2454/' '23X22/' ],
% [experiment_group1 '081618/' '2454/' '23X22/' ],
% [experiment_group1 '082018/' '2454/' '23X22/' ],
% [experiment_group1 '082318/' '2454/' '23X22/' ],
% [experiment_group1 '082718/' '2454/' '23X22/' ],
% [experiment_group1 '083018/' '2454/' '23X22/' ]};
%
% {[experiment_group1 '080618/' '2467/' '23X22/' ],
% [experiment_group1 '080918/' '2467/' '23X22/' ],
% [experiment_group1 '081318/' '2467/' '23X22/' ],
% [experiment_group1 '081618/' '2467/' '23X22/' ],
% [experiment_group1 '082018/' '2467/' '23X22/' ],
% [experiment_group1 '082318/' '2467/' '23X22/' ],
% [experiment_group1 '082718/' '2467/' '23X22/' ],
% [experiment_group1 '083018/' '2467/' '23X22/' ]};
%
% {[experiment_group1 '080618/' '2468/' '23X22/' ],
% [experiment_group1 '080918/' '2468/' '23X22/' ],
% [experiment_group1 '081318/' '2468/' '23X22/' ],
% [experiment_group1 '081618/' '2468/' '23X22/' ],
% [experiment_group1 '082018/' '2468/' '23X22/' ],
% [experiment_group1 '082318/' '2468/' '23X22/' ],
% [experiment_group1 '082718/' '2468/' '23X22/' ],
% [experiment_group1 '083018/' '2468/' '23X22/' ]};
%
% {[experiment_group2 '081518/' '2455/' '23X22/' ],
% [experiment_group2 '082118/' '2455/' '23X22/' ],
% [experiment_group2 '082418/' '2455/' '23X22/' ],
% [experiment_group2 '082818/' '2455/' '23X22/' ],
% [experiment_group2 '083118/' '2455/' '23X22/' ],
% [experiment_group2 '090418/' '2455/' '23X22/' ],
% };
%
% {[experiment_group2 '081518/' '2465/' '23X22/' ],
% [experiment_group2 '082118/' '2465/' '23X22/' ],
% [experiment_group2 '082418/' '2465/' '23X22/' ],
% [experiment_group2 '082818/' '2465/' '23X22/' ],
% [experiment_group2 '083118/' '2465/' '23X22/' ],
% [experiment_group2 '090418/' '2465/' '23X22/' ],
% };
%
% {[experiment_group2 '081518/' '2469/' '23X22/' ],
% [experiment_group2 '082118/' '2469/' '23X22/' ],
% [experiment_group2 '082418/' '2469/' '23X22/' ],
% [experiment_group2 '082818/' '2469/' '23X22/' ],
% [experiment_group2 '083118/' '2469/' '23X22/' ],
% [experiment_group2 '090418/' '2469/' '23X22/' ]
% };
%
% {[experiment_group3  '190321/' '3722/' 'dense/'],
% [experiment_group3  '190325/' '3722/' 'dense/'],
% [experiment_group3  '190328/' '3722/' 'dense/'],
% %[experiment_group3  '190401/' '3722/' 'dense/'],
% [experiment_group3  '190404/' '3722/' 'fullscan/'],
% [experiment_group3  '190409/' '3722/' 'fullscan/'],
% [experiment_group3  '190415/' '3722/' 'fullscan/']};
%
% {[experiment_group3  '190321/' '3724/' 'dense/'],
% [experiment_group3  '190325/' '3724/' 'dense/'],
% [experiment_group3  '190328/' '3724/' 'dense/'],
% %[experiment_group3  '190401/' '3724/' 'dense/'],
% [experiment_group3  '190404/' '3724/' 'fullscan/'],
% [experiment_group3  '190408/' '3724/' 'fullscan/']
% [experiment_group3  '190415/' '3724/' 'fullscan/']};
%
% {[experiment_group3  '190321/' '3726/' 'dense/'],
% [experiment_group3  '190325/' '3726/' 'dense/'],
% [experiment_group3  '190328/' '3726/' 'dense/'],
% %[experiment_group3  '190401/' '3726/' 'dense/'],
% [experiment_group3  '190404/' '3726/' 'fullscan/'],
% [experiment_group3  '190408/' '3726/' 'fullscan/'],
% [experiment_group3  '190415/' '3726/' 'fullscan/']};
%
% {[experiment_group4  '031219/' '2469/' '23X22/'],
% [experiment_group4  '190318/' '2469/' '23X22/'],
% [experiment_group4  '190324/' '2469/' 'dense/'],
% [experiment_group4  '032719/' '2469/' 'dense/'],
% [experiment_group4  '190406/' '2469/' 'Activity_Scan/'],
% [experiment_group4  '190410/' '2469/' 'dense/'],
% [experiment_group4  '190413/' '2469/' 'Activity_Scan/']};
%
% {[experiment_group5  '190323/' '3723/' 'dense/'],
% [experiment_group5  '190326/' '3723/' 'dense/'],
% [experiment_group5  '190329/' '3723/' 'dense/'],
% [experiment_group5  '190402/' '3723/' 'dense/'],
% [experiment_group5  '190405/' '3723/' 'dense/'],
% [experiment_group5  '190412/' '3723/' 'dense/'],
% [experiment_group5  '190416/' '3723/' 'dense/']};
%
% {[experiment_group5  '190323/' '3725/' 'dense/'],
% [experiment_group5  '190326/' '3725/' 'dense/'],
% [experiment_group5  '190329/' '3725/' 'dense/'],
% [experiment_group5  '190402/' '3725/' 'dense/'],
% %[experiment_group5  '190405/' '3725/' 'dense/'],
% [experiment_group5  '190412/' '3725/' 'dense/'],
% [experiment_group5  '190416/' '3725/' 'dense/']}
% }


%
% rec_types = {
%     {'2squares', '2squares', '2squares','2squares', '2squares', '2squares', '2squares', '2squares'};
%     {'2squares', '2squares', '2squares','2squares', '2squares', '2squares', '2squares', '2squares'};
%     {'2squares', '2squares', '2squares','2squares', '2squares', '2squares', '2squares', '2squares'};
%     {'2squares','2squares', '2squares','2squares', '2squares', '2squares'};
%     {'2squares','2squares', '2squares','2squares', '2squares', '2squares'};
%     {'2squares','2squares', '2squares','2squares', '2squares', '2squares'};
%     {'2squares', '2squares', '2squares','1squares', '1squares', '1squares'};
%     {'2squares', '2squares', '2squares','1squares', '1squares', '1squares'};
%     {'2squares', '2squares', '2squares','1squares', '1squares', '1squares'};
%     {'2squares', '2squares', '2squares','2squares','1squares','2squares','1squares'};
%     {'2squares','2squares', '2squares','2squares', '2squares', '2squares', '2squares'  };
%     {'2squares','2squares', '2squares','2squares', '2squares', '2squares' };
%     }





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% analysis options %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pipeline options
spikeDetectionCompleted = true;
corrCalcCompleted = false;
syncCalcCompleted = false;

n_workers = 4;

%spike detection options
n_chunks = 6 * 10;
chunk_size = round(10 * 20000);
list_of_thresholds = [8];
selected_threshold_idx = 1; %index of threshold to use in further calculations
selected_threshold = list_of_thresholds(selected_threshold_idx);

%covariance calculation options
binning_factors = [1, 5].*20;
order_factor = 10;
realMaxOrders = order_factor .* binning_factors;
CMatricesFolderNameBase = 'CMatrices_binned_';
pValuesFolderNameBase = 'pValues_binned_binary_';


%synchronized regions options
activity_threshold = 0.1 %Hz %0.01
selected_binning_factors = [1, 5].*20; %determine synchronicity for these binning factors (in ms)
pvalue_alpha = 0.01; %used in the statistical test to determine significant pairs
sync_thresholds = [0, 0.25, 0.75]; %pairs with zero lag cross-correlations above this threshold are considered synchronized


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% analysis of dense recordings %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for date_idx = 1:length(dates)
    
    selected_date = dates{date_idx};
    %parfor (chip_idx = 1:length(chips), n_workers)
    for chip_idx = 1:length(chips)
        
        
        %        %parfor (d_idx = 1:length(data_to_read), n_workers)
        % for d_idx = 1:length(data_to_read)
        %     for p_idx = 1:length(data_to_read{d_idx})
        %
        %                 resultpath = [global_path data_to_read{d_idx}{p_idx}]
        %                 basepath = [global_path basepaths{d_idx}];
        %
        %
        %                 rec_info = importdata([ resultpath 'rec.csv']);
        %                 sec_in_rec = rec_info(:,2);
        %                 min_spikes = activity_threshold*sec_in_rec;
        %
        %
        %                 files = dir(resultpath);
        %                 dirFlags = [files.isdir];
        %                 subdirs = files(dirFlags);
        %                 recording_resultpaths = {};
        %                 for m = 1:length(subdirs)
        %                     if ~isempty(regexp(subdirs(m).name,'[0-9]+','match'))
        %                         recording_resultpaths{end+1} = [resultpath subdirs(m).name];
        %                     end
        %                 end
        %
        %                 squares = [1];
        %
        %                 if strcmp(rec_types{d_idx}{p_idx}, '2squares')
        %                     squares = [1, 2];
        %                 elseif strcmp(rec_types{d_idx}{p_idx}, '1squares')
        %                     squares = [1];
        %                 end
        %
        %
        %                 path_parts = strsplit(resultpath,filesep)
        %                 date = path_parts(8)
        %                 date = date{1};
        %                 chip = path_parts(9)
        %                 chip = chip{1};
        
        
        
        
        %%%%%%%%% prepare for processing %%%%%%%%%%
        
        chip = chips{chip_idx};
        array_of_matching_chips_files = array_of_matching_files(contains(array_of_matching_files, ['/' chip '/']))
        
        list_of_recordings = {};
        recording_resultpaths = {};
        
        
        for k = 1:length(array_of_matching_chips_files)
            
            filename = array_of_matching_chips_files(k)
            filename = filename{1};
            
            path_parts = strsplit(filename,filesep);
            date =  path_parts(6);
            date = date{1}
            %date = selected_date;
            
            recnum_and_extension_splitted = split(path_parts(end),'.');
            rec = str2num(recnum_and_extension_splitted{1});
            
            if strcmp(date, selected_date)
                
                recording_file = filename
                basepath = [global_path];
                resultpath = [basepath selected_date '/' chip '/' rec_type '/'];
                mkdir(resultpath)
                %basepath = [global_path experiment_group3];
                %resultpath = [basepath date '/' chip '/' 'fullscan' '/'];
                
                
                %%%%% some preprocessing %%%%%%
                
                
                
                %                    list_of_recordings = {};
                %                    recording_resultpaths = {};
                n_channels_per_rec = [];
                n_channels_cumulative = [];
                total_n_channels = 0;
                
                %for rec_idx = 1:length(recordings)
                %for rec_idx = 1:length(recording_resultpaths)
                %rec = recordings(rec_idx);
                %recording_file = [datapath basename num2str(rec, recNumberFormat) filextension]
                %recording_file = [datapath basename filextension]
                
                %recording_file = recording_resultpaths{rec_idx}
                
                list_of_recordings{end+1} = recording_file;
                
                recording_resultpath = [resultpath num2str(rec, recNumberFormat) '/'];
                mkdir(recording_resultpath);
                recording_resultpaths{end+1} = recording_resultpath;
                
                recording = mxw.fileManager(recording_file);
                [testdata, ~, ~] = recording.extractRawData(1, 10);
                total_n_channels = total_n_channels + size(testdata,2);
                n_channels_per_rec = [n_channels_per_rec size(testdata,2)];
                n_channels_cumulative = [n_channels_cumulative total_n_channels];
                
                %%%%%% plot electrode locations %%%%
                
                map = recording.rawMap.map;
                for square = squares
                    
                    if (length(squares) == 1) % single square recording, or network recording
                        square_x = 1:length(map.x); %all indices in map
                        
                    else %recording consisting of two squares
                        rel_middle = (max(map.x) - min(map.x))/2;
                        if (square == 1)
                            square_x = map.x < (min(map.x) + rel_middle); %indices in map correspoding to first square
                        else
                            square_x = map.x >= (min(map.x) + rel_middle);%indices in map correspoding to second square
                        end
                    end

                    f = figure();
                    scatter(map.x(square_x), map.y(square_x), 10);
                    ylim([0 2500]);
                    xlim([0 4000]);
                    pictureName = ['electrodes_locations_' int2str(square)];
                    path = [recording_resultpath '/' int2str(square) '/'];
                    mkdir(path);
                    savepng( 'Directory', path, 'FileName' , pictureName );
                    hold off;
                    close(f)
                    dlmwrite([path 'positions.csv' ],[map.x(square_x), map.y(square_x)] );
                    
                end     
            end
        end
        
       
        
        
        %%%%%%%%% spike detection %%%%%%%%%%
        
        if ~spikeDetectionCompleted
            
            spiketimes_per_recording = {};
            amplitudes_per_recording = {};
            
            total_nspikes = []; %list of n_spikes detected for all sensors (for total summary statistics)
            total_mean_amplitudes = [];  %list of mean amplitude detected for all sensors (for total summary statistics)
            
            
            total_x = [];
            total_y = [];
            
            
            for r = 1:length(list_of_recordings)
                
                disp('started detecting spikes in recording...')
                recording_name = list_of_recordings(r)
                resultpath = [basepath selected_date '/' chip '/' rec_type '/'];
                recording_resultpath = recording_resultpaths(r);
                
                recording = mxw.fileManager(recording_name{1});
                map = recording.rawMap.map;
                total_x = [total_x transpose(map.x)];
                total_y = [total_y transpose(map.y)];
                
                %[all_threshold_spiketimes,all_threshold_amplitudes] = detect_spikes(list_of_recordings, chunk_size, n_chunks, list_of_thresholds, resultpath, recording_resultpaths);
                [threshold_spiketimes,threshold_amplitudes] = detect_spikes(recording_name, chunk_size, n_chunks, list_of_thresholds, resultpath, recording_resultpath);
                
                
                % for now, only use the first given threshold
                threshold_idx = 1;
                
                spiketimes = threshold_spiketimes{threshold_idx}; %spiketimes_per_recording{end};
                amplitudes = threshold_amplitudes{threshold_idx}; %amplitudes_per_recording{end};
                
                
                spiketimes_per_recording{end+1} = spiketimes;
                amplitudes_per_recording{end+1} = amplitudes;
                
                
                for sensor_idx = 1:size(spiketimes,2)
                    
                    sensor_spiketimes = spiketimes(1,sensor_idx);
                    sensor_spiketimes = sensor_spiketimes{1};
                    total_nspikes = [total_nspikes length(sensor_spiketimes)];
                    
                    if length(sensor_spiketimes) >= 1
                        
                        sensor_amplitudes = amplitudes(1,sensor_idx);
                        sensor_amplitudes = sensor_amplitudes{1};
                        total_mean_amplitudes = [total_mean_amplitudes  mean(sensor_amplitudes)];
                    else
                        total_mean_amplitudes = [total_mean_amplitudes  0];
                    end
                end
                
            end
            
            %%% plot and save summary data %%%%
            
            disp('plotting summaries...')
            
            rec_info = importdata([resultpath 'rec.csv']);
            sec_in_rec = rec_info(:,2);
            
            % for now, only use the first given threshold
            threshold = list_of_thresholds(1);
            
            figure()
            histogram(total_nspikes/sec_in_rec);
            title(['Detected Spike Frequency, with threshold ' int2str(threshold)])
            xlabel('Frequency (Hz)')
            ylabel('Number of Sensors')
            pictureName = ['sensor_spikes_histogram_' int2str(threshold)];
            savepng( 'Directory', resultpath, 'FileName' , pictureName );
            hold off
            
            figure()
            histogram(abs(total_mean_amplitudes));
            title(['Detected mean spike amplitude, with threshold ' int2str(threshold)])
            xlabel('Abs. Value Amplitude (microV)')
            ylabel('Number of Sensors')
            pictureName = ['sensor_amplitudes_histogram_' int2str(threshold)];
            savepng( 'Directory', resultpath, 'FileName' , pictureName );
            hold off
            
            
            nspikes_csv = [resultpath 'sensor_nspikes_'  int2str(threshold)  '.csv'];
            meanAmpl_csv = [resultpath 'sensor_meanAmpl_'  int2str(threshold)  '.csv'];
            %since we append to this file, check whether it should be
            %removed -> so we don't append to an older version of the file
            if exist(nspikes_csv, 'file') == 2
                system(['rm -rf ' nspikes_csv]);
            end
            
            if exist(meanAmpl_csv, 'file') == 2
                system(['rm -rf ' meanAmpl_csv]);
            end
            
            
            for k = 1:length(total_nspikes)
                dlmwrite(nspikes_csv, total_nspikes(k), '-append','delimiter',',');
                
                if(total_nspikes(k)~=0)
                    dlmwrite(meanAmpl_csv, total_mean_amplitudes(k), '-append','delimiter',',');
                    %else
                    %   dlmwrite(meanAmpl_csv, 0, '-append','delimiter',',');
                end
                
            end
            
            
            
            
            
            
            %             %split total spiketimes structure in spiketimes per recording
            %             spiketimes_per_recording{end+1} = all_threshold_spiketimes{selected_threshold_idx}(1:n_channels_cumulative(1,1));
            %             amplitudes_per_recording{end+1} = all_threshold_amplitudes{selected_threshold_idx}(1:n_channels_cumulative(1,1));
            %
            %             for j = 2:length(n_channels_per_rec)
            %                 start = n_channels_cumulative(1,j-1)+1;%n_channels_cumulative(1,j-1);
            %                 stop = n_channels_cumulative(1,j);
            %                 spiketimes_per_recording{end+1} = all_threshold_spiketimes{selected_threshold_idx}(start:stop);
            %                 amplitudes_per_recording{end+1} = all_threshold_amplitudes{selected_threshold_idx}(start:stop);
            %             end
            
            
            %%%%%%%%%%%%%% activity maps %%%%%%%%%%%%%%%%%
            
            for t = 1:length(list_of_thresholds)
                
                threshold = list_of_thresholds(t);
                
                nspikes = importdata([resultpath 'sensor_nspikes_' int2str(threshold) '.csv']);
                mean_ampl = importdata([resultpath 'sensor_meanAmpl_' int2str(threshold) '.csv']);
                rec_info = importdata([resultpath 'rec.csv']);
                %                 meta = importdata([resultpath 'meta.csv']);
                %                 x_list =  meta(:,2);
                %                 y_list = meta(:,3);
                sec_in_rec = rec_info(:,2);
                
                cmapLoad = load('cmap_bluered.mat')
                mycmap = cmapLoad.mycmap;
                f = figure();
                colormap(mycmap./256);
                clims = [-2,2];
                plot_2D_map_clean(transpose(total_x), transpose(total_y), total_nspikes/sec_in_rec, clims, 'nearest');
                title(['Mean Detected Spike Rate per Sensor, with threshold ' int2str(threshold) ' MAD']);
                xlabel('\mu m ');
                ylabel('\mu m ');
                colorbar;
                outpath =  [resultpath '/activity_maps_freq/'];
                mkdir(outpath)
                pictureName = [date '_' chip '_freq_' int2str(threshold)];
                savepng( 'Directory', outpath, 'FileName' , pictureName );
                hold off;
                close(f);
                
                
                f=figure();
                colormap(mycmap./256);
                clims = [-100,100];
                plot_2D_map_clean(transpose(total_x), transpose(total_y), total_mean_amplitudes, clims, 'nearest');
                title(['Mean Detected Amplitude per Sensor, with threshold ' int2str(threshold) ' MAD'])
                xlabel('\mu m ');
                ylabel('\mu m ');
                colorbar;
                outpath =  [resultpath '/activity_maps_ampl/'];
                mkdir(outpath)
                pictureName = [date '_' chip '_ampl_'  int2str(threshold)];
                savepng( 'Directory', outpath, 'FileName' , pictureName );
                hold off;
                close(f);
                
            end
            
            
            %%%%%%%%%%% analysis and plots per square %%%%%%%%%%%%%
            
            for i = 1:length(spiketimes_per_recording)
                
                recording = mxw.fileManager(list_of_recordings{i});
                map = recording.rawMap.map;
                
                for square = squares
                    
                    if (length(squares) == 1) % single square recording, or network recording
                        square_idx = 1:length(map.x); %all indices in map
                    else %recording consisting of two squares
                        rel_middle = (max(map.x) - min(map.x))/2;
                        if (square == 1)
                            square_idx = map.x < (min(map.x) + rel_middle); %indices in map corresponding to first square
                        else
                            square_idx = map.x >= (min(map.x) + rel_middle);%indices in map corresponding to second square
                        end
                    end
                    
                    %save the indices in the original data per square
                    dlmwrite([recording_resultpaths{i} '/'  int2str(square) '/' 'indices_in_data.csv' ],square_idx);
                    
                    %save the spiketimes per square in a .mat file
                    sensor_spiketimes = spiketimes_per_recording{i}(square_idx);
                    filenameSpiketimes = ['spiketimes_' int2str(selected_threshold) '.mat'];
                    resultfile_mat = [ [recording_resultpaths{i} '/'  int2str(square) '/'] filenameSpiketimes ];
                    parsave(resultfile_mat,sensor_spiketimes)
                    
                    %save the amplitudes per square in a .mat file
                    sensor_amplitudes = amplitudes_per_recording{i}(square_idx);
                    filenameAmplitudes = ['amplitudes_' int2str(selected_threshold) '.mat'];
                    resultfile_mat = [ [recording_resultpaths{i} '/'  int2str(square) '/']  filenameAmplitudes];
                    parsave(resultfile_mat,sensor_amplitudes)
                    
                    
                    %create analysis plots
                    binning_factor_activity_plot = 1000; %number of samples
                    ymax = 500;
                    binned_activity_plot(spiketimes_per_recording{i}(square_idx), binning_factor_activity_plot, selected_threshold, 20000, ymax, [recording_resultpaths{i} '/'  int2str(square) '/']);
                    %binned_activity_plot(amplitudes_per_recording{i}(square_idx), binning_factor_activity_plot, selected_threshold, 20000, ymax, [recording_resultpaths{i} '/'  int2str(square) '/']);
                    rasterPlotwithAmplitudes(spiketimes_per_recording{i}(square_idx), amplitudes_per_recording{i}(square_idx),  [recording_resultpaths{i} '/'  int2str(square) '/'], selected_threshold);
                    
                end
            end
        end
        
        
        rec_info = importdata([ resultpath 'rec.csv']);
        sec_in_rec = rec_info(:,2)
        min_spikes = activity_threshold*sec_in_rec
        
        nsecond_filedata = importdata([resultpath 'rec.csv']);
        n_seconds_processed = nsecond_filedata(1,2);
        
        
        %calculate covariances
        if ~corrCalcCompleted
            disp('starting correlation calculations')
            
            for b =1:length(binning_factors)
                binning_factor = binning_factors(b);
                realMaxOrder = realMaxOrders(b);
                
                for j = 1:length(recording_resultpaths)
                    recording_resultpaths(j)
                    for square = squares
                        squarePath =  [recording_resultpaths{j} '/'  int2str(square) '/'];
                        filenameSpiketimes = ['spiketimes_' int2str(selected_threshold) '.mat'];
                        
                        spiketimeCorrelations(squarePath, filenameSpiketimes, binning_factor, realMaxOrder, n_seconds_processed, CMatricesFolderNameBase,pValuesFolderNameBase);
                        
                        if length(squares) > 1 %for recordings in two squares, where each recording thus covers different electrodes, we can compute the p-values right away
                            spikeTimePValues(squarePath, binning_factor, realMaxOrder,CMatricesFolderNameBase, pValuesFolderNameBase, n_seconds_processed);
                        end
                    end
                end
                
                if length(squares) == 1%the recordings were multiple recordings of fixed electrodes spread across the array. We can add the spike counts (auto and cross) for all the electrodes that were recorded during the different sessions
                    
                    %determine the fixed electrodes and their
                    %indices in the covariance matrices
                    idx_per_rec = {};
                    if rec_type == 'axon_recording' %we can us the total spike count for all the electrodes that were shared across all the recordings
                        
                  
                        %determine the minimal set of fixed_electrodes
                        %shared accross all recordings
                        recording_0 = mxw.fileManager(list_of_recordings{1});
                        map_0 = recording_0.rawMap.map;
                        recording_1 = mxw.fileManager(list_of_recordings{2});
                        map_1 = recording_1.rawMap.map;
                        [fixed_electrodes, ~, ~] = intersect(map_0.electrode, map_1.electrode);
                        
                        for l = 3:length(list_of_recordings)
                            recording_l = mxw.fileManager(list_of_recordings{l});
                            map_l = recording_l.rawMap.map;
                            [~, ~, idxFixed] = intersect(map_l.electrode, fixed_electrodes);
                            fixed_electrodes = fixed_electrodes(idxFixed); 
                            
                        end
                        
                        %determine the indices of these electrodes in each
                        %recording
                                                
                        for l = 1:length(list_of_recordings)
                            recording_l = mxw.fileManager(list_of_recordings{l});
                            map_l = recording_l.rawMap.map;
                            [~, idxl, ~] = intersect(map_l.electrode, fixed_electrodes);
                            idx_per_rec{end+1} = idxl;
                        end
                        
                        recording_1 = mxw.fileManager(list_of_recordings{1});
                        map_1 = recording_1.rawMap.map;
                        fixed_electrode_positions = [map_1.x(idx_per_rec{1}),map_1.y(idx_per_rec{1})];

                        
                    else %each recording contains the same fixed electrodes, and thus all of these electrodes can be used
                        
                        for l = 1:length(list_of_recordings)
                            recording_l = mxw.fileManager(list_of_recordings{l});
                            map_l = recording_l.rawMap.map;
                            idx_per_rec{end+1} = 1:length(map_l.electrode);
                        end
                    end
                    
                    
                    %sum all the spike counts of the fixed electrodes over all the recordings in a
                    %single cov matrix
                    maxOrder = realMaxOrder/binning_factor;
                    total_C_matrices = zeros([length(fixed_electrodes),length(fixed_electrodes),maxOrder]);
                    
                    
                    for j = 1:length(recording_resultpaths)
                        for square = squares
                            squarePath =  [recording_resultpaths{j} '/'  int2str(square) '/'];
                            
                            specific_part_folder_name = [int2str(binning_factor/20) 'ms_' int2str(realMaxOrder/20) 'ms/'];
                            CMatricesFolder = [squarePath CMatricesFolderNameBase specific_part_folder_name];
                            CMatrices = readCMatrices(CMatricesFolder, maxOrder);
                           
                            
                            n_var = size(CMatrices,1);
                            offsets = 0:maxOrder-1;
                            offsets = offsets .* n_var; 
                            
                            for k = 1:maxOrder
                                offset = offsets(k);
                                total_C_matrices(:,:,k) = total_C_matrices(:,:,k) + CMatrices(idx_per_rec{j},idx_per_rec{j}+offset);
                            end
                           
                            
                        end
                    end
                    
          
                    
                    %write this result to a separate folder and add the
                    %pvalue result to this folder as well
                    %the positions of the shared electrodes also need to be stored
                    %here 
                    fixed_electrodes_path = [resultpath '/fixed_electrodes/'];
                    mkdir(fixed_electrodes_path)
                    
                    CMatricesFolder = [fixed_electrodes_path CMatricesFolderNameBase specific_part_folder_name];
                    if ~(exist(CMatricesFolder ,'dir')==7)
                        mkdir(CMatricesFolder)
                    else
                        system(['rm -rf ' CMatricesFolder])
                        mkdir(CMatricesFolder)
                    end
                    
                    for order = 1:maxOrder
                        filename = strcat(CMatricesFolder, 'C_', int2str(order-1), '.csv');
                        dlmwrite(filename,total_C_matrices(:,:,order),'delimiter',',','precision',7);
                        filename = strcat(CMatricesFolder, 'C_', int2str(order-1), '.csv');
                        dlmwrite(filename,total_C_matrices(:,:,order),'delimiter',',','precision',7);
                    end
               
                    spikeTimePValues(fixed_electrodes_path, binning_factor, realMaxOrder,CMatricesFolderNameBase, pValuesFolderNameBase, n_seconds_processed);
                    
                    dlmwrite([fixed_electrodes_path 'positions.csv' ],fixed_electrode_positions );
    
                end
            end
            disp('correlation calculations completed')
        end
        
        
        %caculate sync regions over all recordings
        if ~syncCalcCompleted
            
            %disp('starting sync calculations')
            
            list_of_square_folders = {};
            list_of_square_numbers = {};
            %list_of_datafiles = {};
            list_of_res_folders = {};
            list_of_min_spikes = {};
            
            
            %t = datetime(2019,7,19,15,0,0,0)
            %firstpath =  dir([recording_resultpaths{1} '/' int2str(1) '/rasterplot_colored_2s.png'])
            %firstpath = dir([resultpath '1ms_10ms_sync_25_synchronized_regions_nSensors.csv']);
            if true %firstpath.date < t
                
                if length(squares) > 1
                    for j = 1:length(recording_resultpaths)
                        for square = squares
                            list_of_min_spikes{end+1} = min_spikes;
                            list_of_res_folders{end+1} = recording_resultpaths{j};
                            list_of_square_folders{end+1} = [recording_resultpaths{j} '/' int2str(square) '/'];
                            list_of_square_numbers{end+1} = square;
                            %rec = recordings(j);
                            %recording_file = [datapath basename num2str(rec, recNumberFormat) filextension];
                            %recording_file = [datapath basename  filextension];
                            %list_of_datafiles{end+1} = recording_file;
                        end
                    end
                else
                    list_of_min_spikes{end+1} = min_spikes;
                    list_of_res_folders{end+1} = [resultpath '/fixed_electrodes/'];
                    list_of_square_folders{end+1} = [resultpath '/fixed_electrodes/'];
                    list_of_square_numbers{end+1} = 1;
                end
                
                for selected_binning_factor = selected_binning_factors
                    for sync_threshold = sync_thresholds
                        
                        specific_part_folder_name = [int2str(selected_binning_factor/20) 'ms_' int2str(selected_binning_factor/20 * order_factor) 'ms'];
                        specific_part_filenames = [int2str(selected_binning_factor/20) 'ms_' int2str(selected_binning_factor/20 * order_factor) 'ms_sync_' int2str(sync_threshold*100)];
                        CMatrixFolderName = [CMatricesFolderNameBase specific_part_folder_name '/']
                        pValueFolderName = [pValuesFolderNameBase specific_part_folder_name '/'];
                        
                        synchronized_regions(basepath, date, chip, resultpath, list_of_square_folders,list_of_square_numbers,list_of_res_folders, list_of_min_spikes, selected_threshold, pvalue_alpha, sync_threshold, CMatrixFolderName, pValueFolderName, specific_part_filenames, order_factor);
                        
                    end
                end
                %else
                %disp([recording_resultpaths{1} ' already processed'])
            end
            
        end
        
        
        %                 catch exception
        %                     %chip
        %                     %date
        %                     %processing_failed{end+1} = [chip, date]
        %                     exception
        %                 end
    end
end






% disp('processing failed for ...')
% for i = 1:length(processing_failed)
%     processing_failed{i}
% end

%
%             rel_middle = (max(map.x) - min(map.x))/2;
%             first_square = map.x < (min(map.x) + rel_middle); %indices in map correspoding to first square
%             second_square = map.x >= (min(map.x) + rel_middle);%indices in map correspoding to second square
%
%
%             %save the indices in the original data per square
%             first_square_indices_original = find(map.x < (min(map.x) + rel_middle));
%             dlmwrite([recording_resultpaths{i} '/1/' 'indices_in_data.csv' ],first_square_indices_original);
%             second_square_indices_original = find(map.x >=  (min(map.x) + rel_middle));
%             dlmwrite([recording_resultpaths{i} '/2/' 'indices_in_data.csv' ],second_square_indices_original);
%
%             %save the spiketimes per square in a .mat file
%             sensor_spiketimes = spiketimes_per_recording{i}(first_square);
%             filename = ['spiketimes_' int2str(selected_threshold) '.mat'];
%             resultfile_mat = [ [recording_resultpaths{i} '/1/'] filename];
%             save(resultfile_mat,'sensor_spiketimes')
%
%             sensor_spiketimes = spiketimes_per_recording{i}(second_square);
%             filename = ['spiketimes_' int2str(selected_threshold) '.mat'];
%             resultfile_mat = [ [recording_resultpaths{i} '/2/'] filename];
%             save(resultfile_mat,'sensor_spiketimes')
%
%
%             %save the amplitudes per square in a .mat file
%             sensor_amplitudes = amplitudes_per_recording{i}(first_square);
%             filename = ['amplitudes_' int2str(selected_threshold) '.mat'];
%             resultfile_mat = [ [recording_resultpaths{i} '/1/'] filename];
%             save(resultfile_mat,'sensor_amplitudes')
%
%             sensor_amplitudes = amplitudes_per_recording{i}(second_square);
%             filename = ['amplitudes_' int2str(selected_threshold) '.mat'];
%             resultfile_mat = [ [recording_resultpaths{i} '/2/'] filename];
%             save(resultfile_mat,'sensor_amplitudes')

%             %create analysis plots
%
%             binned_activity_plot(spiketimes_per_recording{i}(first_square), binning_factor, selected_threshold, 20000, ymax, [recording_resultpaths{i} '/1/']);
%             binned_activity_plot(spiketimes_per_recording{i}(second_square), binning_factor, selected_threshold, 20000, ymax, [recording_resultpaths{i} '/2/']);
%
%             rasterPlotwithAmplitudes(spiketimes_per_recording{i}(first_square), amplitudes_per_recording{i}(first_square),  [recording_resultpaths{i} '/1/'], selected_threshold);
%             rasterPlotwithAmplitudes(spiketimes_per_recording{i}(second_square), amplitudes_per_recording{i}(second_square),  [recording_resultpaths{i} '/2/'], selected_threshold);
%
%             %N = 50;
%             %analyseISI(spiketimes_per_recording{i}(first_square), N, binning_factor,  selected_threshold, [recording_resultpaths{i} '/1/']);
%             %analyseISI(spiketimes_per_recording{i}(second_square), N, binning_factor,  selected_threshold, [recording_resultpaths{i} '/2/']);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% analysis of network recordings %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% type= 'network'
% basename = 'network_session0'
%
%
% for date_idx = 1:length(dates)
%
%
%     date = dates{date_idx}
%     datapath = ['/home/euler/MAXONE/dataMaria/' group '/' date '/' chip '/' type '/']
%     resultpath =  ['/home/euler/MAXONE/resultsMaria/' group '/' date '/' chip '/' type '/']
%
%     mkdir(resultpath)
%
%     list_of_recordings = {};
%     recording_file = [datapath basename filextension]
%     list_of_recordings{end+1} = recording_file;
%     recording = mxw.fileManager(recording_file);
%
%     %%%%%% plot electrode locations %%%%
%
%     map = recording.rawMap.map;
%     f=figure()
%     scatter(map.x, map.y, 10)
%     pictureName = ['electrodes_locations'];
%     savepng( 'Directory', resultpath, 'FileName' , pictureName );
%     hold off;
%     close(f)
%
%
%     %%%%%%%%% spike detection %%%%%%%%%%
%     n_chunks = 2 * 6;
%     chunk_size = round(10 * 20000);
%     list_of_thresholds = [8]
%     [all_threshold_spiketimes,all_threshold_amplitudes] = detect_spikes(list_of_recordings, chunk_size, n_chunks, list_of_thresholds, resultpath)
%     selected_threshold_idx = 1;
%     selected_threshold = list_of_thresholds(selected_threshold_idx);
%
%     spiketimes_per_recording = {};
%     amplitudes_per_recording = {};
%     spiketimes_per_recording{end+1} = all_threshold_spiketimes{selected_threshold_idx}; %only a single recording file
%     amplitudes_per_recording{end+1} = all_threshold_amplitudes{selected_threshold_idx};
%     %%%%%%%%%%% plots %%%%%%%%%%
%     binning_factor = 1000; %number of samples
%     ymax = 400;
%     binned_activity_plot(spiketimes_per_recording{1}, binning_factor,selected_threshold, 20000, ymax, resultpath);
%     rasterPlotwithAmplitudes(spiketimes_per_recording{1}, amplitudes_per_recording{1}, resultpath, selected_threshold);
% end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% analysis of random/fixed recordings %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% type= 'axonScan'
% basename = 'axonScan_'
%
% for date_idx = 1:length(dates)
%
%
%     date = dates{date_idx}
%     datapath = ['/home/euler/MAXONE/dataMaria/group1_2454_2467_2468/' date '/' chip '/' type '/']
%     resultpath =  ['/home/euler/MAXONE/resultsMaria/group1_2454_2467_2468/' date '/' chip '/' type '/']
%
%     mkdir(resultpath)
%
%     recordings = 0:9%strsplit(int2str(0:5))
%
%     list_of_recordings = {};
%     recording_resultpaths = {};
%     n_channels_per_rec = [];
%     n_channels_cumulative = [];
%     total_n_channels = 0;
%
%     for rec_idx = 1:length(recordings)
%         rec = recordings(rec_idx)
%         recording_file = [datapath basename sprintf('%03.f',rec) filextension]
%         recording_resultpath = [resultpath sprintf('%03.f',rec) '/'];
%         mkdir(recording_resultpath)
%         list_of_recordings{end+1} =  recording_file;
%         recording_resultpaths{end+1} = recording_resultpath;
%         recording = mxw.fileManager(recording_file);
%         [testdata, ~, ~] = recording.extractRawData(1, 10);
%         total_n_channels = total_n_channels + size(testdata,2);
%         n_channels_per_rec = [n_channels_per_rec size(testdata,2)];
%         n_channels_cumulative = [n_channels_cumulative total_n_channels];
%
%         %%%%%% plot electrode locations %%%%
%
%         map = recording.rawMap.map;
%         figure()
%         scatter(map.x, map.y, 10)
%         pictureName = ['electrodes_locations'];
%         savepng( 'Directory', recording_resultpath, 'FileName' , pictureName );
%         hold off;
%
%     end
%
%
%
%
%     %%%%%%%%% spike detection %%%%%%%%%%
%     n_chunks = 4;
%     chunk_size = round(10 * 20000);
%     list_of_thresholds = [5, 7]
%     [all_threshold_spiketimes,all_threshold_amplitudes] = detect_spikes(list_of_recordings, chunk_size, n_chunks, list_of_thresholds, resultpath)
%
%
%     selected_threshold_idx = 2;
%     selected_threshold = list_of_thresholds(selected_threshold_idx);
%
%
%     spiketimes_per_recording = {};
%
%
%     %split total spiketimes structure in spiketimes per recording
%     all_threshold_spiketimes{selected_threshold_idx}
%     n_channels_cumulative(1,1)
%     all_threshold_spiketimes{selected_threshold_idx}(1:n_channels_cumulative(1,1))
%
%
%     spiketimes_per_recording{end+1} = all_threshold_spiketimes{selected_threshold_idx}(1:n_channels_cumulative(1,1))
%     for j = 2:length(n_channels_per_rec)
%         start = n_channels_cumulative(1,j-1)
%         stop = n_channels_cumulative(1,j)
%         spiketimes_per_recording{end+1} = all_threshold_spiketimes{selected_threshold_idx}(start:stop)
%     end
%
%     %%%%%%%% binned activity plot %%%%%%%
%
%     binning_factor = 1000; %number of samples
%     for i = 1:length(spiketimes_per_recording)
%         ymax = 400;
%         binned_activity_plot(spiketimes_per_recording{i}, binning_factor, selected_threshold, 20000, ymax, recording_resultpaths{i});
%     end
%
%
%
%
%
%
%
%
%
%
% end





%%%%%%%%%%%%%%%%%%%%%%%%%%%

% list_of_thresholds = [8]
%
% for i = 1:length(list_of_thresholds)
%     threshold = list_of_thresholds(i)
%     type = 'network' %'23X22' %'axonScan'
%
%     chips = {'2454', '2467', '2468'}
%     group = 'group1_2454_2467_2468'
%     dates = {'080618', '080918', '081318', '081618', '082018', '082318', '082718', '083018'}
%     timepoints = [11 14 18 21 25 28 32 35]%plating date 07/26/18
%
%
%     mean_amplitudes = zeros(length(chips), length(dates));
%     mean_amplitudes_active = zeros(length(chips), length(dates));
%     mean_nspikes =  zeros(length(chips), length(dates));
%     mean_nspikes_active =  zeros(length(chips), length(dates));
%     nseconds = zeros(length(chips), length(dates));
%
%     for c_idx = 1:length(chips)
%         for d_idx =  1:length(dates)
%
%             chip = chips{c_idx}
%             date = dates{d_idx}
%             resultpath =  ['/home/euler/MAXONE/resultsMaria/' group '/' date '/' chip '/' type '/']
%
%             nspikes = importdata([resultpath 'sensor_nspikes_' int2str(threshold) '.csv']);
%             amplitudes = importdata([resultpath 'sensor_meanAmpl_' int2str(threshold) '.csv']);
%
%             mean_nspikes(c_idx,d_idx) = mean(nspikes);
%             mean_amplitudes(c_idx,d_idx) = abs(mean(amplitudes));
%
%             indices_active = nspikes>0;
%             mean_nspikes_active(c_idx,d_idx) = mean(nspikes(indices_active));
%             mean_amplitudes_active(c_idx,d_idx) = abs(mean(amplitudes(indices_active)));
%
%             nsecond_filedata = importdata([resultpath 'rec.csv'])
%             nseconds(c_idx,d_idx) = nsecond_filedata(1,2);
%
%
%         end
%     end
%
%
%
%
%
%
%
%     chips2 = {'2469', '2465','2455' }
%     group2 = 'group2_2455_2465_2469'
%     dates2 = {'081518', '082118', '082418','082818','083118','090418'};
%     timepoints2 = [9 15 18 22 25 29] %plating date 08/06/18
%
%
%     mean_amplitudes2 = zeros(length(chips2), length(dates2));
%     mean_amplitudes_active2 = zeros(length(chips2), length(dates2));
%     mean_nspikes2 =  zeros(length(chips2), length(dates2));
%     mean_nspikes_active2 =  zeros(length(chips2), length(dates2));
%     nseconds2 = zeros(length(chips2), length(dates2));
%
%     for c_idx = 1:length(chips2)
%         for d_idx =  1:length(dates2)
%
%             chip = chips2{c_idx}
%             date = dates2{d_idx}
%             resultpath =  ['/home/euler/MAXONE/resultsMaria/' group2 '/' date '/' chip '/' type '/']
%
%             nspikes = importdata([resultpath 'sensor_nspikes_' int2str(threshold) '.csv']);
%             amplitudes = importdata([resultpath 'sensor_meanAmpl_' int2str(threshold) '.csv']);
%
%             mean_nspikes2(c_idx,d_idx) = mean(nspikes);
%             mean_amplitudes2(c_idx,d_idx) = abs(mean(amplitudes));
%
%             indices_active = nspikes>0;
%             mean_nspikes_active2(c_idx,d_idx) = mean(nspikes(indices_active));
%             mean_amplitudes_active2(c_idx,d_idx) = abs(mean(amplitudes(indices_active)));
%
%             nsecond_filedata = importdata([resultpath 'rec.csv'])
%             nseconds2(c_idx,d_idx) = nsecond_filedata(1,2);
%
%
%         end
%     end

%
%     resultpath =  ['/home/euler/MAXONE/resultsMaria/' group '/resultsFeb22/']
%     mkdir(resultpath);
%
%
%
%     f = figure()
%     plot(timepoints, mean_nspikes(1,:)./nseconds(1,:), 'yo-')
%     hold on
%     plot(timepoints, mean_nspikes(2,:)./nseconds(2,:), 'ro-')
%     hold on
%     plot(timepoints, mean_nspikes(3,:)./nseconds(3,:), 'bo-')
%     hold on
%     plot(timepoints2, mean_nspikes2(1,:)./nseconds2(1,:), 'yx-')
%     hold on
%     plot(timepoints2, mean_nspikes2(2,:)./nseconds2(2,:), 'rx-')
%     hold on
%     plot(timepoints2, mean_nspikes2(3,:)./nseconds2(3,:), 'bx-')
%
%
%     title(['Mean Detected Spike Rate over 1000 selected sensors, with threshold ' int2str(threshold) ' MAD'])
%     xlabel('Days after Replating ')
%     ylabel('Spike Rate (Hz)')
%     legend('duplication_1', 'control_1', 'deletion_1', 'duplication_2', 'control_2', 'deletion_2')
%     ax = gca;
%     ax.Box = 'off';
%     ax.LineWidth = 2;
%     ax.FontSize = 10;
%     xlim([8 35])
%     xticks(0:35)
%     legend('Location','northwest')
%
%     pictureName = ['sensor_network_freq_' int2str(threshold)];
%     savepng( 'Directory', resultpath, 'FileName' , pictureName );
%     hold off
%     close(f)
%
%
%     f = figure()
%     plot(timepoints, mean_nspikes_active(1,:)./nseconds(1,:), 'o-')
%     hold on
%     plot(timepoints, mean_nspikes_active(2,:)./nseconds(2,:), 'o-')
%     hold on
%     plot(timepoints, mean_nspikes_active(3,:)./nseconds(3,:), 'o-')
%     title(['Mean Detected Spike Rate over selected, active sensors, with threshold ' int2str(threshold) ' MAD'])
%     xlabel('Days after Replating ')
%     ylabel('Spike Rate (Hz)')
%     legend('duplication', 'control', 'deletion')
%     ax = gca;
%     ax.Box = 'off';
%     ax.LineWidth = 2;
%     ax.FontSize = 10;
%     xlim([8 35])
%     xticks(0:35)
%
%     pictureName = ['sensor_network_active_freq_' int2str(threshold)];
%     savepng( 'Directory', resultpath, 'FileName' , pictureName );
%     hold off
%     close(f)
%
%
%
%
%     f = figure()
%     plot(timepoints, mean_amplitudes(1,:), 'o-')
%     hold on
%     plot(timepoints, mean_amplitudes(2,:), 'o-')
%     hold on
%     plot(timepoints, mean_amplitudes(3,:), 'o-')
%     title(['Mean Detected Amplitude over selected sensors, with threshold ' int2str(threshold) ' MAD'])
%     xlabel('Days after Replating ')
%     ylabel('Amplitude (microV)')
%     legend('duplication', 'control', 'deletion')
%     ax = gca;
%     ax.Box = 'off';
%     ax.LineWidth = 2;
%     ax.FontSize = 10;
%     xlim([8 35])
%     xticks(0:35)
%
%     pictureName = ['sensor_network_ampl_' int2str(threshold)];
%     savepng( 'Directory', resultpath, 'FileName' , pictureName );
%     hold off
%     close(f)
%
%
%     f = figure()
%     plot(timepoints, mean_amplitudes_active(1,:), 'o-')
%     hold on
%     plot(timepoints, mean_amplitudes_active(2,:), 'o-')
%     hold on
%     plot(timepoints, mean_amplitudes_active(3,:), 'o-')
%     title(['Mean Detected Amplitude over selected,active sensors, with threshold ' int2str(threshold) ' MAD'])
%     xlabel('Days after Replating ')
%     ylabel('Amplitude (microV)')
%     legend('duplication', 'control', 'deletion')
%     ax = gca;
%     ax.Box = 'off';
%     ax.LineWidth = 2;
%     ax.FontSize = 10;
%     xlim([8 35])
%     xticks(0:35)
%
%     pictureName = ['sensor_network_active_ampl_' int2str(threshold)];
%     savepng( 'Directory', resultpath, 'FileName' , pictureName );
%     hold off
%     close(f)
%
%
%
% end









quit






%%%%%%%%% network bursting %%%%%%%%%