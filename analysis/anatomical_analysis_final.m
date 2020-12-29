%% anatomical_analysis.m
% John Bernabei
% June 2020
% Center for Neuroengineering & Therapeutics

%% set up workspace
clear all

% load required data
load nodeIDsForJohn

load all_regions_raw

all_WM_data =readtable('/Users/jbernabei/Downloads/Atlas_labels/JHU48.txt');
all_WM_data = all_WM_data{:,2};
all_WM_data{49} = 'Temporal WM';
all_WM_data{50} = 'Frontal WM';

all_GM_data = {'Visual';
            'Somatomotor';
            'Dorsal Attention';
            'Ventral Attention';
            'Limbic';
            'Frontoparietal';
            'Default'};

% session list
session_list = {'HUP069','HUP133','HUP133','HUP136','HUP139','HUP140',...
                'HUP142','HUP143','HUP145','HUP146','HUP150','HUP152',...
                'HUP153','HUP154','HUP157','HUP160','HUP165','HUP168',...
                'HUP171','HUP178','HUP179','HUP179','HUP181','HUP182',...
                'HUP187','HUP191','HUP191'};
            
ant_CR = [23,24];
SLF_SFOF = [41,42,43,44];
thal_rad = [29,30];
limbic_wm = [6,35,36,37,38,39,40,45,46];

composite_regions(1).data = ant_CR;
composite_regions(2).data = SLF_SFOF;
composite_regions(3).data = thal_rad;
composite_regions(4).data = limbic_wm;
composite_regions(5).data = 49;
composite_regions(6).data = 50;

% Region list
temporal_list = [1:16,30:31,82:83];
frontal_list = [20:21, 24:29, 50:59, 68:73, 76:81];
occipital_list = [64:67, 22:23];
parietal_list = [32:33, 60:63];
%% Data processing figure for all patients

for polarity = 1:2
    
    if polarity==1
        load graphChLbls
    else 
        load powChLbls
    end

% do basic localization of all nodes
for session = 1:27
    
        % get patient
        patient = session_list{session};
        
        % clear out variables
        clear this_T1_regions
        clear type_ind
        clear elec_mni_coords
        
        % load the electrodes for that session
        
        % Load electrode names & coordinates in MNI space from CSV file
        elec_label_mni_raw = readtable(sprintf('/Users/jbernabei/Documents/PhD_Research/Neurosurgery_SEEG/improved_localization/new_%s_localization.csv',patient));
        
        if polarity==1
            
            ch_names = gchlbl{session};
           
            v = 1;
            for e = 1:length(ch_names)
                elec1 = ch_names{e,1};
                
                
                coord_ind_1 = find(strcmp(elec_label_mni_raw{:,1},elec1));
                
                
                mni_1 = elec_label_mni_raw{coord_ind_1,3:5};
                
                try elec_mni_coords(v,:) = mni_1;
                    seg_roi{v} = elec_label_mni_raw{coord_ind_1,2};
                    v = v+1;
                catch ME
                    fprintf('could not find electrode\n')
                end
            end
        else
            % need to find centroids of bipolar pair
            ch_names = gchlbl{session};
            num_pairs = size(ch_names,1);
            v = 1;
            for e = 1:num_pairs
                elec1 = ch_names{e,1};
                elec2 = ch_names{e,2};
                
                coord_ind_1 = find(strcmp(elec_label_mni_raw{:,1},elec1));
                coord_ind_2 = find(strcmp(elec_label_mni_raw{:,1},elec2));
                
                mni_1 = elec_label_mni_raw{coord_ind_1,3:5};
                mni_2 = elec_label_mni_raw{coord_ind_2,3:5};
                
                try elec_mni_coords(v,:) = (mni_1+mni_2)./2;
                    seg_roi1{v} = elec_label_mni_raw{coord_ind_1,2};
                    seg_roi2{v} = elec_label_mni_raw{coord_ind_2,2};
                    v = v+1;
                catch ME
                    fprintf('could not find electrode\n')
                end
            end
        end
        
        % get number of electrodes
        num_elecs = size(elec_mni_coords,1);
        
        % determine electrode type: 0 = gm, 1 = wm, 2 = out of brain
        type_ind = zeros(num_elecs,1);
        
        if polarity==1
            % Correct ROI to go from blank to 'n/a'
            for e = 1:num_elecs
                if strcmp(seg_roi{e},'')
                    this_T1_regions{e,1} = ('n/a');
                    type_ind(e) = 2;
                else
                    this_T1_regions{e,1} = seg_roi{e};
                end

                if strcmp(this_T1_regions{e,1},'Left Cerebral White Matter')||strcmp(this_T1_regions{e,1},'Right Cerebral White Matter')
                    type_ind(e) = 1;
                end
            end
        else
            % Correct ROI to go from blank to 'n/a'
            for e = 1:num_elecs
                if strcmp(seg_roi1{e},'') && strcmp(seg_roi2{e},'')
                    this_T1_regions{e,1} = ('n/a');
                    type_ind(e) = 2;
                else
                    this_T1_regions{e,1} = seg_roi1{e};
                end

                if strcmp(seg_roi1{e},'Left Cerebral White Matter')||strcmp(seg_roi1{e},'Right Cerebral White Matter')...
                        && strcmp(seg_roi2{e},'Left Cerebral White Matter')||strcmp(seg_roi2{e},'Right Cerebral White Matter')
                    type_ind(e) = 1;
                end
            end
        end
        
        % calculate number of GM & number of WM
        num_gm(session) = sum(type_ind==0);
        num_wm(session) = sum(type_ind==1);
        num_ex_brain(session) = sum(type_ind==2);
        
        % localize gm in gm atlas and wm in wm atlas
        [mni_coords_gm, mni_labels_gm, NN_flag_gm] = nifti_values(elec_mni_coords((type_ind==0),:),'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii');
        [mni_coords_wm, mni_labels_wm, NN_flag_wm] = nifti_values(elec_mni_coords((type_ind==1),:),'/Users/jbernabei/Downloads/Atlas_labels/JHU48.nii');
    
        wm_coords = elec_mni_coords(type_ind==1,:);
        mni_labels_gen = [];
        [mni_coords_gen, mni_labels_gen, NN_flag_gen] = nifti_values(wm_coords,'/Users/jbernabei/Downloads/Hammers_mith_atlas_n30r83_delivery_Dec16/Hammers_mith_atlas_n30r83_SPM5.nii');

        for j = 1:length(mni_labels_wm)
            if mni_labels_wm(j)==0 && ismember(mni_labels_gen(j),temporal_list)
                mni_labels_wm(j)=49;
            elseif mni_labels_wm(j)==0 && ismember(mni_labels_gen(j),frontal_list)
                mni_labels_wm(j)=50;
            end
        end
                        
        num_localized_gm(session) = sum(mni_labels_gm~=0);
        num_localized_wm(session) = sum(mni_labels_wm~=0);
        
%         % put the rest of WM in a general atlas (lobar)
%         wm_coords = elec_mni_coords((type_ind==1),:);
%         mni_labels_gen = [];
%         [mni_coords_gen, mni_labels_gen, NN_flag_gen] = nifti_values(wm_coords(mni_labels_wm==0,:),'/Users/jbernabei/Downloads/Hammers_mith_atlas_n30r83_delivery_Dec16/Hammers_mith_atlas_n30r83_SPM5.nii');
%       
%         for j = 1:length(mni_labels_gen)
%             if ismember(mni_labels_gen(j),temporal_list)
%                 mni_labels_wm(j)=49;
%             elseif ismember(mni_labels_gen(j),frontal_list)
%                 mni_labels_wm(j)=50;
%             end
%         end
        
%         general_wm(session).data = mni_labels_gen;
        
        % gm distribution
        for r = 1:7
            gm_distribution(r,session) = sum(mni_labels_gm==r);
        end
        
        % wm distribution
        for r = 1:50
            wm_distribution(r,session) = sum(mni_labels_wm==r);
        end
        
        % put localization data into this structure
        patient_loc(session).polarity(polarity).labels = elec_label_mni_raw{:,1};
        patient_loc(session).polarity(polarity).coords = elec_mni_coords;
        patient_loc(session).polarity(polarity).T1_loc = this_T1_regions;
        patient_loc(session).polarity(polarity).type = type_ind;
        patient_loc(session).polarity(polarity).yeo = mni_labels_gm;
        patient_loc(session).polarity(polarity).jhu = mni_labels_wm;
        
        for r = 1:7
            region_count_GM(r,session) = sum(patient_loc(session).polarity(polarity).yeo==r);
        end
        
        for r = 1:50
            region_count_WM(r,session) = sum(patient_loc(session).polarity(polarity).jhu==r);
        end
        
        %general_wm(session).data
end


% write the tables 
region_table_GM = table(all_GM_data,      region_count_GM(:,1),  region_count_GM(:,2),  region_count_GM(:,3),...
                 region_count_GM(:,4),   region_count_GM(:,5),  region_count_GM(:,6),  region_count_GM(:,7),...
                 region_count_GM(:,8),   region_count_GM(:,9),  region_count_GM(:,10), region_count_GM(:,11),...
                 region_count_GM(:,12),  region_count_GM(:,13), region_count_GM(:,14), region_count_GM(:,15),...
                 region_count_GM(:,16),  region_count_GM(:,17), region_count_GM(:,18), region_count_GM(:,19),...
                 region_count_GM(:,20),  region_count_GM(:,21), region_count_GM(:,22), region_count_GM(:,23),...
                 region_count_GM(:,24),  region_count_GM(:,25), region_count_GM(:,26), region_count_GM(:,27)); 

             if polarity==1
                %writetable(region_table_GM,'all_gm_counts_monopolar.csv'); 
             else
                % writetable(region_table_GM,'all_gm_counts_bipolar.csv'); 
             end

region_table_WM = table(all_WM_data, region_count_WM(:,1),  region_count_WM(:,2),  region_count_WM(:,3),...
                 region_count_WM(:,4),   region_count_WM(:,5),  region_count_WM(:,6),  region_count_WM(:,7),...
                 region_count_WM(:,8),   region_count_WM(:,9),  region_count_WM(:,10), region_count_WM(:,11),...
                 region_count_WM(:,12),  region_count_WM(:,13), region_count_WM(:,14), region_count_WM(:,15),...
                 region_count_WM(:,16),  region_count_WM(:,17), region_count_WM(:,18), region_count_WM(:,19),...
                 region_count_WM(:,20),  region_count_WM(:,21), region_count_WM(:,22), region_count_WM(:,23),...
                 region_count_WM(:,24),  region_count_WM(:,25), region_count_WM(:,26), region_count_WM(:,27)); 

             
             if polarity==1
                %writetable(region_table_WM,'all_wm_counts_monopolar.csv'); 
             else
                 %writetable(region_table_WM,'all_wm_counts_bipolar.csv'); 
             end


end

%% Do anatomical analysis

qexp_vs_pow = {'qexp','pow'};
threshold = [0.05, 0.1, 1];
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];
load color_bar

qq = 0;

% do for qexp and pow
for metric = 1:2
    
    if metric==1
        load graphChLbls
    else 
        load powChLbls
    end
    
    % do for frequency bands 1:4
    for freq = 1:4
        
        qq = qq+1;
        
        % do for thresholds 0.05 and 0.1
        for thres = 3
            
            for r = 1:13
            all_sig_b_val(r).data = [];
            end
            
            cmap_GM = [];
            cmap_WM = [];
            
            % do for sessions 1:27
            for session = 1:27
                %session
                clear ch_status
                clear ch_coords
                clear mni_labels_gm
                clear mni_labels_wm
                clear sig_ch
                clear mni_labels_max_gen
                clear mni_labels_min_gen
                clear mni_labels_max_GM
                clear mni_labels_min_GM
                clear mni_labels_max_WM
                clear mni_labels_min_WM
                
                % Get all b and p values
                if strcmp(qexp_vs_pow{metric},'qexp')
                    all_b_val = LTqexp{session}(:,freq,1);
                    all_p_val = LTqexp{session}(:,freq,4);
                else
                    all_b_val = LTpow{session}(:,freq,1);
                    all_p_val = LTpow{session}(:,freq,4);
                end
                
                ch_names = gchlbl{session};

                % get significant channels based on p value threshold
                sig_ch = find(all_p_val<threshold(thres));
                
                % check if there are any significant channels
                if isempty(sig_ch)
                    fprintf('skipping session b/c no significance\n')
                    
                else
                    % extract significant channel p values, p values, and names
                    sig_p_val = all_p_val(sig_ch);
                    sig_b_val = all_b_val(sig_ch);
                    sig_ch_names = ch_names(sig_ch,1);
                    
                    [ind_max] =  (sig_b_val==max(sig_b_val));
                    [ind_min] =  (sig_b_val==min(sig_b_val));
                    
                    
                    
                    % match up significant channels using their names to
                    % coordinates and ROI. Eliminate channels that do not have
                    % corresponding ROI
                    for ch = 1:length(sig_ch_names)
                        % get this individual channel name
                        this_ch = sig_ch_names{ch};
                        
                        % find where in the original channel list for
                        % localization this significant channel is
                        coord_index = find(strcmp(patient_loc(session).polarity(metric).labels,this_ch));
                        try patient_coord_indices(ch) = coord_index;
                            ch_status(ch) = patient_loc(session).polarity(metric).type(coord_index);
                            
                            ch_coords(ch,:) = patient_loc(session).polarity(metric).coords(coord_index,:);
                            sig_max_ind(ch) = ind_max(coord_index);
                            sig_min_ind(ch) = ind_min(coord_index);
                        catch ME
                            patient_coord_indices(ch) = 0;
                            ch_status(ch) = 3;
                            ch_coords(ch,:) = [0 0 0];
                        end
                        
                    end
                    
                    [ch_status' ind_max ind_min]
                    
                    % localize max and min
                    try[mni_coords_max, mni_labels_max_GM, NN_flag] = nifti_values(ch_coords(find(sig_max_ind),:),'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii');
                    catch ME
                        mni_labels_max_GM = [];
                    end
                    try [mni_coords_max, mni_labels_max_WM, NN_flag] = nifti_values(ch_coords(find(sig_max_ind),:),'JHU48.nii');
                    catch ME
                        mni_labels_max_WM = [];
                    end
                    try [mni_coords_max, mni_labels_max_gen, NN_flag] = nifti_values(ch_coords(find(sig_max_ind),:),'/Users/jbernabei/Downloads/Hammers_mith_atlas_n30r83_delivery_Dec16/Hammers_mith_atlas_n30r83_SPM5.nii');
                    catch ME
                        mni_labels_max_gen = [];
                    end
                    if mni_labels_max_GM~=0
                        max_region{session,qq} = all_GM_data{mni_labels_max_GM};
                    elseif mni_labels_max_WM~=0
                        max_region{session,qq} = all_WM_data{mni_labels_max_WM};
                    elseif ismember(mni_labels_max_gen,temporal_list)
                        max_region{session,qq} = 'Temporal WM';
                    elseif ismember(mni_labels_max_gen,frontal_list)
                        max_region{session,qq} = 'Frontal WM';
                    else
                        max_region{session,qq} = 'Unlocalized';
                    end
                    
                    % min
                    try [mni_coords_max, mni_labels_min_GM, NN_flag] = nifti_values(ch_coords(find(sig_min_ind),:),'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii');
                    catch ME
                        mni_labels_min_GM = [];
                    end
                    try [mni_coords_max, mni_labels_min_WM, NN_flag] = nifti_values(ch_coords(find(sig_min_ind),:),'JHU48.nii');
                    catch ME
                        mni_labels_min_WM = [];
                    end
                    try [mni_coords_max, mni_labels_min_gen, NN_flag] = nifti_values(ch_coords(find(sig_min_ind),:),'/Users/jbernabei/Downloads/Hammers_mith_atlas_n30r83_delivery_Dec16/Hammers_mith_atlas_n30r83_SPM5.nii');
                    catch ME
                        mni_labels_min_gen = [];
                    end
                    if mni_labels_min_GM~=0
                        min_region{session,qq} = all_GM_data{mni_labels_min_GM};
                    elseif mni_labels_min_WM~=0
                        min_region{session,qq} = all_WM_data{mni_labels_min_WM};
                    elseif ismember(mni_labels_min_gen,temporal_list)
                        min_region{session,qq} = 'Temporal WM';
                    elseif ismember(mni_labels_min_gen,frontal_list)
                        min_region{session,qq} = 'Frontal WM';
                    else
                        min_region{session,qq} = 'Unlocalized';
                    end
                    
                    % extract b, p values and do GM and WM locs on the significant nodes
                    b_val_gm = sig_b_val(find((ch_status==0)+(ch_status==2)));
                    b_val_wm = sig_b_val(find(ch_status==1));
                    
                    mni_labels_gm = [];
                    try [mni_coords_gm, mni_labels_gm, NN_flag_gm] = nifti_values(ch_coords(find(ch_status==0),:),'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii');
                    catch ME
                        fprintf('mni GM problem\n')
                    end
                    
                    % find number of significant nodes in each of the 7
                    % systems
                    for r = 1:7
                        region_count_GM(r,session) = sum(mni_labels_gm==r);
                    end
                    
                    mni_labels_wm = [];
                    try [mni_coords_wm, mni_labels_wm, NN_flag_wm] = nifti_values(ch_coords(ch_status==1,:),'/Users/jbernabei/Downloads/Atlas_labels/JHU48.nii');
                         % put the rest of WM in a general atlas (lobar)
                        wm_coords = ch_coords(ch_status==1,:);
                        mni_labels_gen = [];
                        [mni_coords_gen, mni_labels_gen, NN_flag_gen] = nifti_values(wm_coords,'/Users/jbernabei/Downloads/Hammers_mith_atlas_n30r83_delivery_Dec16/Hammers_mith_atlas_n30r83_SPM5.nii');

                        for j = 1:length(mni_labels_wm)
                            if mni_labels_wm(j)==0 && ismember(mni_labels_gen(j),temporal_list)
                                mni_labels_wm(j)=49;
                            elseif mni_labels_wm(j)==0 && ismember(mni_labels_gen(j),frontal_list)
                                mni_labels_wm(j)=50;
                            end
                        end
                    catch ME
                        fprintf('mni WM problem\n')
                    end
                    
                    % find number of significant nodes in each of the 48
                    % WM regions
                    for r = 1:50
                        region_count_WM(r,session) = sum(mni_labels_wm==r);
                    end
                    
                    % for each of 7 Yeo regions find mean b val in patient
                    for r = 1:7
                        try all_sig_b_val(r).data = [all_sig_b_val(r).data; (b_val_gm(ismember(mni_labels_gm,r)))];
                        catch ME
                            fprintf('nothing in gm region\n')
                            r
                        end
                    end
                    
                    % for each of 4 aggregate WM regions find mean b val in
                    % patient
                    for r = 8:13
                        corrected_r = r-7;
                        try all_sig_b_val(r).data = [all_sig_b_val(r).data; (b_val_wm(ismember(mni_labels_wm,composite_regions(corrected_r).data)))];
                        catch ME
                            fprintf('nothing in wm region\n')
                            corrected_r
                        end
                    end
                   
              
                end
            end
            
%             max_nodes = max([length(all_sig_b_val(1).data),length(all_sig_b_val(2).data),...
%                              length(all_sig_b_val(3).data),length(all_sig_b_val(4).data),...
%                              length(all_sig_b_val(5).data),length(all_sig_b_val(6).data),...
%                              length(all_sig_b_val(7).data),length(all_sig_b_val(8).data),...
%                              length(all_sig_b_val(9).data),length(all_sig_b_val(10).data),...
%                              length(all_sig_b_val(11).data),length(all_sig_b_val(12).data),...
%                              length(all_sig_b_val(13).data)]);
%                          
%             data_table = NaN*ones(max_nodes,13);
%             
%             % fill 11 roi
%             for r = 1:13
%                 num_sig_nodes = length(all_sig_b_val(r).data);
%                 data_table(1:num_sig_nodes,r) = all_sig_b_val(r).data;
%             end
%             
%             data_table_neg = data_table;
%             data_table_pos = data_table;
%             
%             data_table_neg(find(1-(data_table_neg<0))) = NaN;
%             data_table_pos(find(1-(data_table_pos>0))) = NaN;
%             
%             [p,tbl,stats] = kruskalwallis(data_table);
%             c_all = multcompare(stats) 
%             
%             start_sig = c_all(find(c_all(:,6)<0.05),1);
%             stop_sig = c_all(find(c_all(:,6)<0.05),2);
%             
%             % plot 1
%             figure(1);clf;
%             hold on
%             boxplot(data_table,'Labels',{'Visual','SM','D. Attention',...
%                                          'V. Attention','Limbic GM','FP',...
%                                          'Default','Ant. CR','SLF-SFOF','Thal. Rad.',...
%                                          'Limbic WM','Temporal WM','Frontal WM'})
%             set(gcf, 'Position',  [100, 100, 800, 500])
%             ylabel('Regression slope')
%             title(sprintf('b value - %s frequency %d, p = %f',qexp_vs_pow{metric},freq,p));
%             h = findobj(gca,'Tag','Box');
%             colors = [color1; color1;color1; color1; color1; color1; color2; color2;color2; color2;color2; color2;color2];
%             for j=1:length(h)
%                 patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
%             end
%             if ~isempty(start_sig)
%                 num_sig = length(start_sig);
%                 for j = 1:num_sig
%                 y_val = max([data_table(:,start_sig(j));data_table(:,stop_sig(j))]);
%                 plot([start_sig(j),stop_sig(j)], [y_val*1.1.^j,y_val*1.1.^j], '-k', 'LineWidth',2)
%                 plot(mean([start_sig(j),stop_sig(j)]), y_val*(1.1.^j+0.05), '*k')
%                 end
%             end
%             hold off
%             %saveas(gcf,sprintf('b_value_%s_frequency_%d.jpg',qexp_vs_pow{metric},freq))
% 
%             % table pos
%             [p,tbl,stats] = kruskalwallis(data_table_pos);
%             c_pos = multcompare(stats) 
%             start_sig = c_pos(find(c_pos(:,6)<0.05),1);
%             stop_sig = c_pos(find(c_pos(:,6)<0.05),2);
%             
%             figure(2);clf;
%             hold on
%             boxplot(data_table_pos,'Labels',{'Visual','SM','D. Attention',...
%                                          'V. Attention','Limbic GM','FP',...
%                                          'Default','Ant. CR','SLF-SFOF','Thal. Rad.',...
%                                          'Limbic WM','Temporal WM','Frontal WM'})
%             set(gcf, 'Position',  [100, 100, 800, 500])
%             ylabel('Regression slope')
%             title(sprintf('positive b value - %s frequency %d, p = %f',qexp_vs_pow{metric},freq,p));
%             h = findobj(gca,'Tag','Box');
%             colors = [color1; color1;color1; color1; color1; color1; color2; color2;color2; color2;color2; color2;color2];
%             for j=1:length(h)
%                 patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
%             end
%             if ~isempty(start_sig)
%                 num_sig = length(start_sig);
%                 for j = 1:num_sig
%                 y_val = max([data_table(:,start_sig(j));data_table(:,stop_sig(j))]);
%                 plot([start_sig(j),stop_sig(j)], [y_val*1.1.^j,y_val*1.1.^j], '-k', 'LineWidth',2)
%                 plot(mean([start_sig(j),stop_sig(j)]), y_val*(1.1.^j+0.05), '*k')
%                 end
%             end
%             hold off
%             %saveas(gcf,sprintf('pos_b_value_%s_frequency_%d.jpg',qexp_vs_pow{metric},freq))
%             
%             % table neg
%             [p,tbl,stats] = kruskalwallis(data_table_neg);
%             c_neg = multcompare(stats) 
%             
%             start_sig = c_neg(find(c_neg(:,6)<0.05),1);
%             stop_sig = c_neg(find(c_neg(:,6)<0.05),2);
%             
%             figure(3);clf;
%             hold on
%             boxplot(data_table_neg,'Labels',{'Visual','SM','D. Attention',...
%                                          'V. Attention','Limbic GM','FP',...
%                                          'Default','Ant. CR','SLF-SFOF','Thal. Rad.',...
%                                          'Limbic WM','Temporal WM','Frontal WM'})
%             set(gcf, 'Position',  [100, 100, 800, 500])
%             ylabel('Regression slope')
%             title(sprintf('negative b value - %s frequency %d, p = %f',qexp_vs_pow{metric},freq,p));
%             h = findobj(gca,'Tag','Box');
%             colors = [color1; color1;color1; color1; color1; color1; color2; color2;color2; color2;color2; color2;color2];
%             for j=1:length(h)
%                 patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
%             end
%             if ~isempty(start_sig)
%                 num_sig = length(start_sig);
%                 for j = 1:num_sig
%                 y_val = max([data_table(:,start_sig(j));data_table(:,stop_sig(j))]);
%                 plot([start_sig(j),stop_sig(j)], [y_val*1.1.^j,y_val*1.1.^j], '-k', 'LineWidth',2)
%                 plot(mean([start_sig(j),stop_sig(j)]), y_val*(1.1.^j+0.05), '*k')
%                 end
%             end
%             hold off
%             %saveas(gcf,sprintf('neg_b_value_%s_frequency_%d.jpg',qexp_vs_pow{metric},freq))
%             
%           
%           
%            % write the tables 
%            region_table_GM = table(all_GM_data,      region_count_GM(:,1),  region_count_GM(:,2),  region_count_GM(:,3),...
%                              region_count_GM(:,4),   region_count_GM(:,5),  region_count_GM(:,6),  region_count_GM(:,7),...
%                              region_count_GM(:,8),   region_count_GM(:,9),  region_count_GM(:,10), region_count_GM(:,11),...
%                              region_count_GM(:,12),  region_count_GM(:,13), region_count_GM(:,14), region_count_GM(:,15),...
%                              region_count_GM(:,16),  region_count_GM(:,17), region_count_GM(:,18), region_count_GM(:,19),...
%                              region_count_GM(:,20),  region_count_GM(:,21), region_count_GM(:,22), region_count_GM(:,23),...
%                              region_count_GM(:,24),  region_count_GM(:,25), region_count_GM(:,26), region_count_GM(:,27)); 
% 
%            %writetable(region_table_GM,sprintf('%s_gm_counts_freq_%d_%f_threshold.csv',qexp_vs_pow{metric},freq,threshold(thres))); 
% 
%            region_table_WM = table(all_WM_data, region_count_WM(:,1),  region_count_WM(:,2),  region_count_WM(:,3),...
%                              region_count_WM(:,4),   region_count_WM(:,5),  region_count_WM(:,6),  region_count_WM(:,7),...
%                              region_count_WM(:,8),   region_count_WM(:,9),  region_count_WM(:,10), region_count_WM(:,11),...
%                              region_count_WM(:,12),  region_count_WM(:,13), region_count_WM(:,14), region_count_WM(:,15),...
%                              region_count_WM(:,16),  region_count_WM(:,17), region_count_WM(:,18), region_count_WM(:,19),...
%                              region_count_WM(:,20),  region_count_WM(:,21), region_count_WM(:,22), region_count_WM(:,23),...
%                              region_count_WM(:,24),  region_count_WM(:,25), region_count_WM(:,26), region_count_WM(:,27)); 

           %writetable(region_table_WM,sprintf('%s_wm_counts_freq_%d_%f_threshold.csv',qexp_vs_pow{metric},freq,threshold(thres))); 

%            [p,tbl,stats] = kruskalwallis([all_gm_b_val(:,1:7),all_wm_b_val]);
%            p_value_metric(metric).threshold(thres).data(freq) = p;
%            
%            c = multcompare(stats);
%            if sum(c(:,6)<0.05)
%                fprintf('significance found\n')
%            end

           % do rendering
%            region_values = nanmedian(data_table);
%            min_val = min(region_values);
%            max_val = max(region_values);
%            
%            if abs(min_val)>abs(max_val)
%                entry = floor(region_values./abs(min_val)*50)+51;
%            else
%                entry = floor(region_values./abs(max_val)*50)+50;
%            end
%            
%            color_map_regions = color_bar(entry,:);
%            
%            cmap_GM = color_map_regions(1:7,:);
%            cmap_WM_raw = color_map_regions(8:11,:);
%            
%            cmap_WM = zeros(17,3);
%            
%            cmap_WM(1,:) = cmap_WM_raw(1,:);
%            cmap_WM(2,:) = cmap_WM_raw(1,:);
%            
%            cmap_WM(3,:) = cmap_WM_raw(2,:);
%            cmap_WM(4,:) = cmap_WM_raw(2,:);
%            cmap_WM(5,:) = cmap_WM_raw(2,:);
%            cmap_WM(6,:) = cmap_WM_raw(2,:);
%            
%            cmap_WM(7,:) = cmap_WM_raw(3,:);
%            cmap_WM(8,:) = cmap_WM_raw(3,:);
%            
%            cmap_WM(9,:) = cmap_WM_raw(4,:);
%            cmap_WM(10,:) = cmap_WM_raw(4,:);
%            cmap_WM(11,:) = cmap_WM_raw(4,:);
%            cmap_WM(12,:) = cmap_WM_raw(4,:);
%            cmap_WM(13,:) = cmap_WM_raw(4,:);
%            cmap_WM(14,:) = cmap_WM_raw(4,:);
%            cmap_WM(15,:) = cmap_WM_raw(4,:);
%            cmap_WM(16,:) = cmap_WM_raw(4,:);
%            cmap_WM(17,:) = cmap_WM_raw(4,:);
%            
%            cmap_WM
%            cmap_GM
%            
           % now we must map colors into composite regions
           %dlmwrite(sprintf('cmap_GM_metric_%f_freq_%f.txt',metric,freq),cmap_GM,'delimiter',' ','precision',5)
           %dlmwrite(sprintf('cmap_WM_metric_%f_freq_%f.txt',metric,freq),cmap_WM,'delimiter',' ','precision',5)
           
           
           %BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','JHU48.nii','JHU_surface_rendering.mat',sprintf('WM_surface_metric_%d_freq_%d.jpg',metric,freq))
           %BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii','Yeo_surface_rendering.mat',sprintf('GM_surface_metric_%d_freq_%d.jpg',metric,freq))
         end
    end
end



%%

for freq = 1:4
    % plot 1: mean b value for different regions
    mean_bval_temp_GM_all = [];
    mean_bval_temp_WM_all = [];
    mean_bval_front_GM_all = [];
    mean_bval_front_WM_all = [];

    
    v = 0;
    b = 0;
    clear sig_roi;
    
    region_count = [];
    region_table = [];
    all_regions = [];
    
    all_sig_b_val_temp_GM = [];
    all_sig_b_val_temp_WM = [];
    all_sig_b_val_front_GM = [];
    all_sig_b_val_front_WM = [];

    
    % loop through sessions
    for session = [1:27]
        session
        
        % get patient from list of sessions
        patient = session_list{session};
            
        % clear out ROIs and coordinates
        clear elec_mni_coords
        clear elec_roi
        clear sig_coords
        clear sig_regions_raw
        clear all_regions_raw
        
        sig_b_val = [];
        sig_p_val = [];

        % Get channel names
        ch_names = gchlbl{session};

        % Get all b and p values
        if strcmp(qexp_vs_pow,'qexp')
            all_b_val = LTqexp{session}(:,freq,1);
            all_p_val = LTqexp{session}(:,freq,4);
        else
            all_b_val = LTpow{session}(:,freq,1);
            all_p_val = LTpow{session}(:,freq,4);
        end
        
        % get significant channel indices and data from those channels
        if strcmp(all_vs_minmax,'all')
            if strcmp(threshold,'low')
                sig_ch = find(all_p_val<0.05);
            else
                sig_ch = find(all_p_val<0.1);
            end
        else
            [~, sig_ch_max] = max(all_b_val);
            [~, sig_ch_min] = min(all_b_val);
            sig_ch = [sig_ch_max, sig_ch_min];
        end
        
        if isempty(sig_ch)
        else
        
        % extract significant channel p values, p values, and names
        sig_p_val = all_p_val(sig_ch);
        sig_b_val = all_b_val(sig_ch);
        sig_ch_names = ch_names(sig_ch,1);
        
        % Load electrode names & coordinates in MNI space from CSV file
        elec_label_mni_raw = readtable(sprintf('/Users/jbernabei/Documents/PhD_Research/Neurosurgery_SEEG/improved_localization/new_%s_localization.csv',patient));
        
        % Get mni coordinates
        elec_mni_coords = elec_label_mni_raw{:,3:5};
        
        clear this_T1_regions
        
        sig_coords = [];
        this_T1_regions = {};
        
        % get coords of sig channels
        for c = 1:length(sig_ch_names)
            this_sig_ch = sig_ch_names{c};
            coord_index = find(strcmp(elec_label_mni_raw{:,1},this_sig_ch));
            
            % get the coordinates we want
            sig_coords(c,1:3) = elec_mni_coords(coord_index,1:3);
            v = v+1;
            if strcmp(elec_label_mni_raw{coord_index,2},'')
                this_T1_regions(c,1) = {'n/a'};
                T1_regions(v,1) = {'n/a'};
            else
                this_T1_regions(c,1) = elec_label_mni_raw{coord_index,2};
                T1_regions(v,1) = elec_label_mni_raw{coord_index,2};
            end
        end
 
        % assemble session data
        session_data = [sig_coords,-1*sig_b_val, 0.1./sig_p_val];
        
        % remove regions that are localized outside of the brain
        remove_regions = find(strcmp(this_T1_regions,'n/a'));
        session_data(remove_regions,:) = [];  
        sig_coords(remove_regions,:) = [];
        sig_b_val(remove_regions) = [];
        sig_p_val(remove_regions) = [];
   
    % correct coords of significant nodes
    sig_coords(find(sig_coords(:,1)>88),1) = 88;
    sig_coords(find(sig_coords(:,2)>107),2) = 107;
    sig_coords(find(sig_coords(:,3)>88),3) = 88;
    sig_coords(find(sig_coords(:,1)<-88),1) = -88;
    sig_coords(find(sig_coords(:,2)<-107),2) = -107;
    sig_coords(find(sig_coords(:,3)<-88),3) = -88;
    
    % correct coords of all nodes
    elec_mni_coords(find(elec_mni_coords(:,1)>88),1) = 88;
    elec_mni_coords(find(elec_mni_coords(:,2)>107),2) = 107;
    elec_mni_coords(find(elec_mni_coords(:,3)>88),3) = 88;
    elec_mni_coords(find(elec_mni_coords(:,1)<-88),1) = -88;
    elec_mni_coords(find(elec_mni_coords(:,2)<-107),2) = -107;
    elec_mni_coords(find(elec_mni_coords(:,3)<-88),3) = -88;
    
    if isempty(sig_coords)
    else
    [mni_coords3, mni_labels3, NN_flag3] = nifti_values(sig_coords,'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii');
    [mni_coords4, mni_labels4, NN_flag4] = nifti_values(elec_mni_coords,'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii');
    [mni_coords5, mni_labels5, NN_flag5] = nifti_values(sig_coords,'/Users/jbernabei/Downloads/Atlas_labels/JHU48.nii');
    [mni_coords6, mni_labels6, NN_flag6] = nifti_values(elec_mni_coords,'/Users/jbernabei/Downloads/Atlas_labels/JHU48.nii');
    
    %a = 1
    
    for r = 1:48
        region_count_WM(r,session) = sum(mni_labels5==r);
        all_regions_WM(r,session) = sum(mni_labels6==r);
    end
    
    for r = 1:7
        region_count_GM(r,session) = sum(mni_labels3==r);
        all_regions_GM(r,session) = sum(mni_labels4==r);
    end
    
        sig_regions_raw = [];
        all_regions_raw = [];
        for i = 1:size(mni_labels3,2)
        if sum((strcmp(this_T1_regions{i,1},'Left Cerebral White Matter') || strcmp(this_T1_regions{i,1},'Right Cerebral White Matter')))
            sig_regions_raw(i) = 1;
        else 
            sig_regions_raw(i) = 2;
        end
        end
        
        for i = 1:length(elec_mni_coords(:,1))
        if sum((strcmp(elec_label_mni_raw{i,2},'Left Cerebral White Matter') || strcmp(elec_label_mni_raw{i,2},'Right Cerebral White Matter')))
            all_regions_raw(i) = 1;
        else 
            all_regions_raw(i) = 2;
        end
        end
        
        num_elec_seg_wm(session,freq) = sum(all_regions_raw==1);
        num_elec_seg_gm(session,freq) = sum(all_regions_raw==2);
        
        num_elec_mni_gm(session,freq) = sum(all_regions_GM(:,session));
        num_elec_mni_wm(session,freq) = sum(all_regions_WM(:,session));
        
        frac_wm_represented(session,freq) = sum(num_elec_mni_wm(session,freq)./num_elec_seg_wm(session,freq))
        frac_gm_represented(session,freq) = sum(num_elec_mni_gm(session,freq)./num_elec_seg_gm(session,freq))
        
        %         
%         
%         % do same for all electrodes in that patient
%         for i = 1:size(mni_labels4,2)
%         if sum(mni_labels4(i)==temporal_list) && (strcmp(elec_label_mni_raw{i,2},'Left Cerebral White Matter') || strcmp(elec_label_mni_raw{i,2},'Right Cerebral White Matter'))
%             all_regions_raw(i) = 2;
%         elseif sum(mni_labels4(i)==frontal_list) && (strcmp(elec_label_mni_raw{i,2},'Left Cerebral White Matter') || strcmp(elec_label_mni_raw{i,2},'Right Cerebral White Matter'))
%             all_regions_raw(i) = 4;
%         elseif sum(strcmp(elec_label_mni_raw{i,2},'Left Cerebral White Matter') || strcmp(elec_label_mni_raw{i,2},'Right Cerebral White Matter'))
%             sig_regions_raw(i) = 5;
%         elseif sum(mni_labels4(i)==temporal_list) 
%             all_regions_raw(i) = 1;
%         elseif sum(mni_labels4(i)==frontal_list) 
%             all_regions_raw(i) = 3;
%         end
%         end
%         
%         
%         temp_GM_inds = find(sig_regions_raw==1);
%         temp_WM_inds = find(sig_regions_raw==2);
%         front_GM_inds = find(sig_regions_raw==3);
%         front_WM_inds = find(sig_regions_raw==4);
%         
%        	non_included_WM_inds_sig(session,freq) = sum(sig_regions_raw==5);
%         included_WM_inds_sig(session,freq) = sum(sig_regions_raw~=5);
%         non_included_WM_inds_all(session) = sum(all_regions_raw==5);
%         included_WM_inds_all(session) = sum(all_regions_raw~=5);
%         
%         mean_bval_temp_GM = nanmean(sig_b_val(temp_GM_inds));
%         mean_bval_temp_WM = nanmean(sig_b_val(temp_WM_inds));
%         mean_bval_front_GM = nanmean(sig_b_val(front_GM_inds));
%         mean_bval_front_WM = nanmean(sig_b_val(front_WM_inds));
% 
%         expected_frac_front_GM = sum(all_regions_raw==3)./(sum(all_regions_raw==1)+sum(all_regions_raw==3));
%         expected_frac_front_WM = sum(all_regions_raw==4)./(sum(all_regions_raw==2)+sum(all_regions_raw==4));
%         
%         all_sig_b_val_temp_GM = [all_sig_b_val_temp_GM; sig_b_val(temp_GM_inds)];
%         all_sig_b_val_temp_WM = [all_sig_b_val_temp_WM; sig_b_val(temp_WM_inds)];
%         all_sig_b_val_front_GM = [all_sig_b_val_front_GM; sig_b_val(front_GM_inds)];
%         all_sig_b_val_front_WM = [all_sig_b_val_front_WM; sig_b_val(front_WM_inds)];
%   
%    % add together stuff from all sessions
%    % figure 1
%     mean_bval_temp_GM_all = [mean_bval_temp_GM_all; mean_bval_temp_GM];
%     mean_bval_temp_WM_all = [mean_bval_temp_WM_all; mean_bval_temp_WM];
%     mean_bval_front_GM_all = [mean_bval_front_GM_all; mean_bval_front_GM];
%     mean_bval_front_WM_all = [mean_bval_front_WM_all; mean_bval_front_WM];
%     
%     observed_frac_front_GM_pos_all = [observed_frac_front_GM_pos_all; observed_frac_front_GM_pos];
%     observed_frac_front_WM_pos_all = [observed_frac_front_WM_pos_all; observed_frac_front_WM_pos];
%     observed_frac_front_GM_neg_all = [observed_frac_front_GM_neg_all; observed_frac_front_GM_neg];
%     observed_frac_front_WM_neg_all = [observed_frac_front_WM_neg_all; observed_frac_front_WM_neg];
%     
        end
        end
    end
   
    region_table_WM = table(all_WM_data{:,2},region_count_WM(:,1),region_count_WM(:,2),region_count_WM(:,3),...
                     region_count_WM(:,4),region_count_WM(:,5),region_count_WM(:,6),region_count_WM(:,7),...
                     region_count_WM(:,8),region_count_WM(:,9),region_count_WM(:,10),region_count_WM(:,11),...
                     region_count_WM(:,12),region_count_WM(:,13),region_count_WM(:,14),region_count_WM(:,15),...
                     region_count_WM(:,16),region_count_WM(:,17),region_count_WM(:,18),region_count_WM(:,19),...
                     region_count_WM(:,20),region_count_WM(:,21),region_count_WM(:,22),region_count_WM(:,23),...
                     region_count_WM(:,24),region_count_WM(:,25),region_count_WM(:,26),region_count_WM(:,27)); 

    writetable(region_table_WM,sprintf('%s_wm_counts_freq_%d_%s_threshold.csv',qexp_vs_pow,freq,threshold)); 
       
    all_region_table_WM = table(all_WM_data{:,2},all_regions_WM(:,1),all_regions_WM(:,2),all_regions_WM(:,3),...
                         all_regions_WM(:,4),all_regions_WM(:,5),all_regions_WM(:,6),all_regions_WM(:,7),...
                         all_regions_WM(:,8),all_regions_WM(:,9),all_regions_WM(:,10),all_regions_WM(:,11),...
                         all_regions_WM(:,12),all_regions_WM(:,13),all_regions_WM(:,14),all_regions_WM(:,15),...
                         all_regions_WM(:,16),all_regions_WM(:,17),all_regions_WM(:,18),all_regions_WM(:,19),...
                         all_regions_WM(:,20),all_regions_WM(:,21),all_regions_WM(:,22),all_regions_WM(:,23),...
                         all_regions_WM(:,24),all_regions_WM(:,25),all_regions_WM(:,26),all_regions_WM(:,27)); 
                     
    writetable(all_region_table_WM,'all_regions_wm.csv'); 
    
    region_table_GM = table(all_GM_data,     region_count_GM(:,1),  region_count_GM(:,2),  region_count_GM(:,3),...
                     region_count_GM(:,4),   region_count_GM(:,5),  region_count_GM(:,6),  region_count_GM(:,7),...
                     region_count_GM(:,8),   region_count_GM(:,9),  region_count_GM(:,10), region_count_GM(:,11),...
                     region_count_GM(:,12),  region_count_GM(:,13), region_count_GM(:,14), region_count_GM(:,15),...
                     region_count_GM(:,16),  region_count_GM(:,17), region_count_GM(:,18), region_count_GM(:,19),...
                     region_count_GM(:,20),  region_count_GM(:,21), region_count_GM(:,22), region_count_GM(:,23),...
                     region_count_GM(:,24),  region_count_GM(:,25), region_count_GM(:,26), region_count_GM(:,27)); 

    writetable(region_table_GM,sprintf('%s_gm_counts_freq_%d_%s_threshold.csv',qexp_vs_pow,freq,threshold)); 
       
    all_region_table_GM = table(all_GM_data,      all_regions_GM(:,1),  all_regions_GM(:,2),  all_regions_GM(:,3),...
                          all_regions_GM(:,4),    all_regions_GM(:,5),  all_regions_GM(:,6),  all_regions_GM(:,7),...
                          all_regions_GM(:,8),    all_regions_GM(:,9),  all_regions_GM(:,10), all_regions_GM(:,11),...
                          all_regions_GM(:,12),   all_regions_GM(:,13), all_regions_GM(:,14), all_regions_GM(:,15),...
                          all_regions_GM(:,16),   all_regions_GM(:,17), all_regions_GM(:,18), all_regions_GM(:,19),...
                          all_regions_GM(:,20),   all_regions_GM(:,21), all_regions_GM(:,22), all_regions_GM(:,23),...
                          all_regions_GM(:,24),   all_regions_GM(:,25), all_regions_GM(:,26), all_regions_GM(:,27)); 
                     
    writetable(all_region_table_GM,'all_regions_gm.csv'); 
   
%     mean_bval_temp_GM_all_freq = [mean_bval_temp_GM_all_freq,mean_bval_temp_GM_all];
%     mean_bval_temp_WM_all_freq = [mean_bval_temp_WM_all_freq,mean_bval_temp_WM_all];
%     mean_bval_front_GM_all_freq = [mean_bval_front_GM_all_freq,mean_bval_front_GM_all];
%     mean_bval_front_WM_all_freq = [mean_bval_front_WM_all_freq,mean_bval_front_WM_all]; 
%     
%     observed_frac_front_GM_pos_all_freq = [observed_frac_front_GM_pos_all_freq, observed_frac_front_GM_pos_all];
%     observed_frac_front_WM_pos_all_freq = [observed_frac_front_WM_pos_all_freq, observed_frac_front_WM_pos_all];
%     observed_frac_front_GM_neg_all_freq = [observed_frac_front_GM_neg_all_freq, observed_frac_front_GM_neg_all];
%     observed_frac_front_WM_neg_all_freq = [observed_frac_front_WM_neg_all_freq, observed_frac_front_WM_neg_all];
%       
        
%         final_plot_GM = [all_sig_mni_GM,-1*all_sig_bval_GM,0.1./all_sig_pval_GM];
%         final_plot_WM = [all_sig_mni_WM,-1*all_sig_bval_WM,0.1./all_sig_pval_WM];
%         
%          dlmwrite('mni_final_render.node',final_plot_GM,'delimiter',' ','precision',5)
%          BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','mni_final_render.node','svm_render.mat',sprintf('%s_render_freq_%d_GM.jpg_pos',qexp_vs_pow,freq))
%           dlmwrite('mni_final_render.node',final_plot_WM,'delimiter',' ','precision',5)
%          BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','mni_final_render.node','svm_render.mat',sprintf('%s_render_freq_%d_WM.jpg_pos',qexp_vs_pow,freq))
%          
%          final_plot_GM_neg = [all_sig_mni_GM_neg,-1*all_sig_bval_GM_neg,0.1./all_sig_pval_GM_neg];
%         final_plot_WM_neg = [all_sig_mni_WM_neg,-1*all_sig_bval_WM_neg,0.1./all_sig_pval_WM_neg];
%         
%          dlmwrite('mni_final_render.node',final_plot_GM_neg,'delimiter',' ','precision',5)
%          BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','mni_final_render.node','svm_render.mat',sprintf('%s_render_freq_%d_GM_neg.jpg',qexp_vs_pow,freq))
%           dlmwrite('mni_final_render.node',final_plot_WM_neg,'delimiter',' ','precision',5)
%          BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','mni_final_render.node','svm_render.mat',sprintf('%s_render_freq_%d_WM_neg.jpg',qexp_vs_pow,freq))

end

%% make plots
for i = 1:4
    p1 = ranksum(observed_frac_front_GM_pos_all_freq(:,i),expected_frac_front_GM);
    p2 = ranksum(observed_frac_front_GM_pos_all_freq(:,i),expected_frac_front_WM);
    p3 = ranksum(observed_frac_front_GM_neg_all_freq(:,i),expected_frac_front_GM);
    p4 = ranksum(observed_frac_front_GM_neg_all_freq(:,i),expected_frac_front_WM);
    
    pvals(:,i) = [p1;p2;p3;p4];
end

%% make plots
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];
color3 = [78 172 91]/255;
color4 = [103 55 155]/255;

for i = 1:4
    plot_data = [mean_bval_temp_GM_all_freq(:,i),...
        mean_bval_temp_WM_all_freq(:,i),...
        mean_bval_front_GM_all_freq(:,i),...
        mean_bval_front_WM_all_freq(:,i)];
    if i==1
        band = 'alpha/theta'
    elseif i==2
        band = 'beta'
    elseif i==3
        band = 'low-gamma'
    else
        band = 'high-gamma'
    end
    
    p1 = ranksum(mean_bval_temp_GM_all_freq(:,i),mean_bval_temp_WM_all_freq(:,i))
    p2 = ranksum(mean_bval_front_GM_all_freq(:,i),mean_bval_front_WM_all_freq(:,i))
    
    p3 = ranksum(mean_bval_temp_GM_all_freq(:,i),mean_bval_front_GM_all_freq(:,i))
    p4 = ranksum(mean_bval_temp_WM_all_freq(:,i),mean_bval_front_WM_all_freq(:,i))
    
    p5 = ranksum(mean_bval_temp_GM_all_freq(:,i),0)
    p6 = ranksum(mean_bval_front_GM_all_freq(:,i),0)
    
    p7 = ranksum(mean_bval_temp_WM_all_freq(:,i),0)
    p8 = ranksum(mean_bval_front_WM_all_freq(:,i),0)
    
    figure(i);clf;
    boxplot(plot_data,'Labels',{'Temporal GM','Temporal WM','Frontal GM','Frontal WM'})
    title(sprintf('Regional regression slope in %s band: Qexp',band))
    colors = [color4; color3; color2; color1];
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
    %legend('Temporal GM','Temporal WM','Frontal GM','Frontal WM','Location','SouthEast')
    
    pvals(:,i) = [p1;p2;p3;p4];

end


%% Make the plots
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];
color3 = [78 172 91]/255;
color4 = [103 55 155]/255;
color5 = [255 179 71]/255;

for i = 1:4
    p1(i) = signrank(Temporal_WM_all(:,i),Temporal_WM_real)
    p2(i) = signrank(Temporal_GM_all(:,i),Temporal_GM_real)
end

for i = 1:4
    p3(i) = signrank(Temporal_WM_all_neg(:,i),Temporal_WM_real_neg)
    p4(i) = signrank(Temporal_GM_all_neg(:,i),Temporal_GM_real_neg)
end

figure(1);clf;
hold on
title('Qexp: Frontal lobe white matter')
boxplot([1-[Temporal_WM_all,NaN*ones(26,1),Temporal_WM_real]],'Labels',{'Alpha/Theta','Beta','Low-Gamma','High-Gamma',' ','Expected'})
ylabel('Fraction of significant nodes in WM')
ylim([-0.1,1.2])
h = findobj(gca,'Tag','Box');
colors = [color5; color1; color4; color3; color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

hold off

figure(2);clf;
hold on
title('Qexp: Frontal lobe grey matter')
boxplot([1-[Temporal_GM_all,NaN*ones(26,1),Temporal_GM_real]],'Labels',{'Alpha/Theta','Beta','Low-Gamma','High-Gamma',' ','Expected'})
ylabel('Fraction of significant nodes in GM')
ylim([-0.1,1.2])
h = findobj(gca,'Tag','Box');
colors = [color5; color1; color4; color3; color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
hold off

figure(3);clf;
hold on
title('Qexp: Frontal lobe white matter')
boxplot([1-[Temporal_WM_all_neg,NaN*ones(26,1),Temporal_WM_real_neg]],'Labels',{'Alpha/Theta','Beta','Low-Gamma','High-Gamma',' ','Expected'})
ylabel('Fraction of significant nodes in WM')
ylim([-0.1,1.2])
h = findobj(gca,'Tag','Box');
colors = [color5; color1; color4; color3; color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end

hold off

figure(4);clf;
hold on
title('Qexp: Frontal lobe grey matter')
boxplot([1-[Temporal_GM_all_neg,NaN*ones(26,1),Temporal_GM_real_neg]],'Labels',{'Alpha/Theta','Beta','Low-Gamma','High-Gamma',' ','Expected'})
ylabel('Fraction of significant nodes in GM')
ylim([-0.1,1.2])
h = findobj(gca,'Tag','Box');
colors = [color5; color1; color4; color3; color2; color1];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
hold off