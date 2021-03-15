%% anatomical_analysis.m
% John Bernabei
% Dec 2020
% Center for Neuroengineering & Therapeutics

%% set up workspace
clear all

% load required data
%load nodeIDsForJohn
load powRT_sIall_10mspre
load graphRT_sIall_10mspre
load CCDTnodeAnalyze

load all_regions_raw

% get all ROI from EVE
EVE_roi_csv = readtable('EVE_all_roi.csv');
EVE_roi_WM_csv = readtable('EVE_regions_WM.csv');
EVE_roi_list = EVE_roi_csv{:,2}; % extract all ROI

% session list - 27 sessions representing 24 patients with SEEG
session_list = {'HUP069','HUP133','HUP133','HUP136','HUP139','HUP140',...
                'HUP142','HUP143','HUP145','HUP146','HUP150','HUP152',...
                'HUP153','HUP154','HUP157','HUP160','HUP165','HUP168',...
                'HUP171','HUP178','HUP179','HUP179','HUP181','HUP182',...
                'HUP187','HUP191','HUP191'};
            
%% quadrant analysis: define quadrants
pow_quad = zeros(size(fPv));
qexp_quad = zeros(size(fQv));

pow_quad(find((fPv>0).*(sPv>0))) = 1;
pow_quad(find((fPv<0).*(sPv>0))) = 2;
pow_quad(find((fPv<0).*(sPv<0))) = 3;
pow_quad(find((fPv>0).*(sPv<0))) = 4;

qexp_quad(find((fQv>0).*(sQv>0))) = 1;
qexp_quad(find((fQv<0).*(sQv>0))) = 2;
qexp_quad(find((fQv<0).*(sQv<0))) = 3;
qexp_quad(find((fQv>0).*(sQv<0))) = 4;
      
%% create EVE roi list bilateral
% we're taking all regions besides 'background' and combining R + L
EVE_roi_list_bilat{1,1} = 'Background';
for i = 2:89
    split_regions = split(EVE_roi_list{i},'_left'); % some string parsing
    EVE_roi_list_bilat(i,1) = split_regions(1);
end
    
% list which ROI of the Hammer smith atlas belong to which lobe (F/T)
temporal_list = [1:16,30:31,82:83];
frontal_list = [20:21, 24:29, 50:59, 68:73, 76:81];

% list top segmentation regions of GM to use based on prevalence
seg_gm_list = {'MFG','MTG','ITG','fusiform','Hippocampus'};

thalamo_cortical = [36,37,38,6,7,35,96,97,98,68,69,95];
frontal_association = [3,4,5,42,43,46,21,24,65,66,67,102,103,106,82,85];
temporal_association = [20,19,18,45,44,12,8,1,81,70,79,105,104,74,70,64];
paralimbic_wm = [41,39,40,11,47,101,100,99,73,107];
commissural = [51,52,53,111,112,113];

frontal_cortex = {'SFG', 'MFG', 'IFG', 'precentral', 'GRe', 'MOrG', 'POrG'};
temporoparietal_cortex = {'MTG', 'ITG', 'STG', 'FuG', 'AnG', 'SMG', 'SPL','postcentral'};
paralimbic_GM = {'Hippocampus','Amygdala','PHG','cingulate','LiG','Ent','insula'};

frontal_cortex_eve = [3,4,5,6,22,24,91,92,93,94,109,112];
temporoparietal_cortex_eve = [1,7,8,12,18,19,20,23,89,95,96,100,106,107,108,111];
paralimbic_gm_eve = [2,11,13,17,25,26,27,90,99,101,105,113,114,115];

% all_region_inds = [frontal_cortex_eve,temporoparietal_cortex_eve,paralimbic_gm_eve]
% [c,i] = sort(all_region_inds)
% all_gm_inds = [ones(1,12),2*ones(1,14),3*ones(1,14)]
% gm_order = all_gm_inds(i)
all_sub_roi(1).data= frontal_cortex_eve;
all_sub_roi(2).data= temporoparietal_cortex_eve;
all_sub_roi(3).data= paralimbic_gm_eve;
all_sub_roi(4).data= thalamo_cortical;
all_sub_roi(5).data= frontal_association;
all_sub_roi(6).data= temporal_association;
all_sub_roi(7).data= paralimbic_wm;
all_sub_roi(8).data= commissural;
%% Electrode processing and localization for all patients

for polarity = 1:2 
    
    % need different labels and thus localization for power vs qexp
    if polarity==1 
        load graphRT_sIall_inst0% monopolar for graph
        polar_type = 'Monopolar';
        quad_nodes_q = [];
        quad_seg_q = [];
        quad_wm_q = [];
    else 
        load powRT_sIall_inst0 % bipolar for power
        polar_type = 'Bipolar';
        quad_nodes_p = [];
        quad_seg_p = [];
        quad_wm_p = [];
    end 
    
    % T1 segmentation ROI structure containing all ROI across patients
    all_T1_roi = [];
 
    % do basic localization of all nodes
    for session = 1:27
    
        % get patient
        patient = session_list{session};
        
        % clear out variables
        clear this_T1_regions
        clear type_ind
        clear elec_mni_coords
        clear seg_roi
        clear seg_roi1
        clear seg_roi2
        clear lobe_list
        clear gm_seg_inds
        clear JHU_wm_list
        
        clear node_type
        clear lobe_type
        clear seg_type
        clear wm_type
        
        % Load electrode names & coordinates in MNI space from CSV file
        elec_label_mni_raw = readtable(sprintf('/Users/jbernabei/Documents/PhD_Research/Neurosurgery_SEEG/improved_localization/new_%s_localization.csv',patient));
        
        % different process for mono vs bipolar montage
        if polarity==1
            
            % extract channel labels
            ch_names = NFstruct(session).gchlbl;
            final_ch_names = ch_names;
            
            v = 1;
            % loop through channels and assign ROI & coords
            for e = 1:length(ch_names)
                
                % extract electrode name
                elec1 = ch_names{e,1};
                
                % match up names
                coord_ind_1 = find(strcmp(elec_label_mni_raw{:,1},elec1));
                
                % extract MNI coordinates
                mni_1 = elec_label_mni_raw{coord_ind_1,3:5};
                
                % get the MNI coords and segmentation ROI
                try elec_mni_coords(v,:) = mni_1;
                    seg_roi(v,1) = elec_label_mni_raw{coord_ind_1,2};
                    v = v+1;
                catch ME
                    % print if electrode is not available 
                    % (segmentation error where a small number of contacts are unlocalized)
                    fprintf('could not find electrode\n')
                    final_ch_names(e,:) = [];
                end
            end
        else
            % however if it is bipolar and not monopolar we 
            % need to find centroids of bipolar pair
            ch_names = NFstruct(session).gchlbl; % get channel names
            final_ch_names = ch_names;
            num_pairs = size(ch_names,1); % find number of bipolar pairs we need
            v = 1;
            
           
            
            % loop through number of channel bipolar pairs
            for e = 1:num_pairs
                
                % get channel names
                elec1 = ch_names{e,1}; 
                elec2 = ch_names{e,2};
                
                
                % get channel indices
                coord_ind_1 = find(strcmp(elec_label_mni_raw{:,1},elec1));
                coord_ind_2 = find(strcmp(elec_label_mni_raw{:,1},elec2));
                
                % get coordinates
                mni_1 = elec_label_mni_raw{coord_ind_1,3:5};
                mni_2 = elec_label_mni_raw{coord_ind_2,3:5};
                
                % assign coordinates as centroid
                try elec_mni_coords(v,:) = (mni_1+mni_2)./2;
                    seg_roi1(v,1) = elec_label_mni_raw{coord_ind_1,2};
                    seg_roi2(v,1) = elec_label_mni_raw{coord_ind_2,2};
                    v = v+1;
                catch ME
                    % print if unavailable
                    fprintf('could not find electrode\n')
                    final_ch_names(e,:) = [];
                end
            end
        end
        
        % get number of electrodes (or pairs in case of bipolar)
        num_elecs = size(elec_mni_coords,1);
        
        % determine electrode type: 0 = gm, 1 = wm, 2 = out of brain
        type_ind = zeros(num_elecs,1);
        
        % need slightly different correction for monopolar vs bipolar
        if polarity==1
            % Correct ROI to go from blank to 'n/a'
            for e = 1:num_elecs
                if strcmp(seg_roi{e},'')
                    this_T1_regions(e,1) = {'n/a'};
                    type_ind(e) = 0;
                else
                    this_T1_regions(e,1) = seg_roi(e);
                    type_ind(e) = 1;
                end

                if strcmp(this_T1_regions{e,1},'Left Cerebral White Matter')||strcmp(this_T1_regions{e,1},'Right Cerebral White Matter')
                    type_ind(e) = -1;
                end
            end
        else
            % Correct ROI to go from blank to 'n/a'
            for e = 1:num_elecs
                if strcmp(seg_roi1{e},'') && strcmp(seg_roi2{e},'')
                    this_T1_regions(e,1) = {'n/a'};
                    type_ind(e) = 0;
                else
                    this_T1_regions(e,1) = seg_roi1(e);
                    type_ind(e) = 1;
                end

                if strcmp(seg_roi1{e},'Left Cerebral White Matter')||strcmp(seg_roi1{e},'Right Cerebral White Matter')...
                        && strcmp(seg_roi2{e},'Left Cerebral White Matter')||strcmp(seg_roi2{e},'Right Cerebral White Matter')
                    type_ind(e) = -1;
                end
            end
        end
        
        gm_seg_inds = zeros(length(this_T1_regions),1);
        
        for r_gm = 1:length(gm_seg_inds)
            this_region = this_T1_regions{r_gm};
            if contains(this_T1_regions{r_gm},frontal_cortex)
                gm_seg_inds(r_gm) = 1;
            elseif contains(this_T1_regions{r_gm},temporoparietal_cortex)
                gm_seg_inds(r_gm) = 2;
            elseif contains(this_T1_regions{r_gm},paralimbic_GM)
                gm_seg_inds(r_gm) = 3;
            end
        end
        
        gm_seg_label_ID = zeros(length(this_T1_regions),1);
        for r_gm = 1:length(gm_seg_label_ID)
            this_region = this_T1_regions{r_gm};
            
            % frontal
            for r_atl = 1:length(frontal_cortex)
                if contains(this_T1_regions{r_gm},frontal_cortex{r_atl})
                    gm_seg_label_ID(r_gm) = r_atl;
                end
            end
            
            % temporoparietal
            for r_atl = 1:length(temporoparietal_cortex)
                if contains(this_T1_regions{r_gm},temporoparietal_cortex{r_atl})
                    gm_seg_label_ID(r_gm) = r_atl+length(frontal_cortex);
                end
            end
            
            % paralimbic
            for r_atl = 1:length(paralimbic_GM)
                if contains(this_T1_regions{r_gm},paralimbic_GM{r_atl})
                    gm_seg_label_ID(r_gm) = r_atl+length(frontal_cortex)+length(temporoparietal_cortex);
                end
            end
        end
        
        % calculate number of GM & number of WM
        num_gm(session) = sum(type_ind==1);
        num_wm(session) = sum(type_ind==-1);
        num_ex_brain(session) = sum(type_ind==0);
        
        % need to render electrodes (type determines color)
%         final_elec_matrix = [elec_mni_coords,type_ind,ones(length(type_ind),1)];
%         dlmwrite('render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
%         BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','render_elecs.node','ecog_seeg_render_2.mat',sprintf('output/Final/Coord_render/%s/elecs_ident_%s.jpg',polar_type,patient))
%         delete render_elecs.node

        % do rendering of WM and GM separately
        % take only GM first
%         gm_coords = elec_mni_coords(find(type_ind==1),:);
%         gm_names = ch_names(find(type_ind==1),:);
        usable_inds = (type_ind==1)+(type_ind==-1);
        
        % JHU EVE atlas - GM+WM
        [mni_coords, node_roi, NN_flag] = nifti_values(elec_mni_coords(usable_inds==1,:),'JHU_MNI_SS_WMPM_Type-I_to_MNI_brain.nii');
        
        % JHU EVE atlas - WM only
        [mni_wm_only, EVE_wm_only, NN_wm_only] = nifti_values(elec_mni_coords(usable_inds==1,:),'JHU_MNI_SS_WMPM_Type-III_to_MNI_brain.nii');
        
        % JHU48 atlas
        [mni_wm_JHU, JHU_wm_only, NN_wm_JHU] = nifti_values(elec_mni_coords(usable_inds==1,:),'JHU48.nii');
        
        % Hammer smith atlas
        [mni_lobar, lobar_roi, NN_lobar] = nifti_values(elec_mni_coords(usable_inds==1,:),'/Users/jbernabei/Downloads/Hammers_mith_atlas_n30r83_delivery_Dec16/Hammers_mith_atlas_n30r83_SPM5.nii');
        
        % extract which lobe each lobe belongs to
        for l = 1:length(lobar_roi)
            if ismember(lobar_roi(l),temporal_list)
                lobe_list(l) = 2;
            elseif ismember(lobar_roi(l),frontal_list)
                lobe_list(l) = 1;
            else
               lobe_list(l) = 0;
            end
        end
        
        %JHU_wm_list = zeros(length(node_wm_only),1);
        for l = 1:length(EVE_wm_only)
            this_wm_roi = EVE_wm_only(l);
            if ismember(this_wm_roi,thalamo_cortical)
                JHU_wm_list(l,1) = 1; 
            elseif ismember(this_wm_roi,frontal_association)
                JHU_wm_list(l,1) = 2;
            elseif ismember(this_wm_roi,temporal_association)
                JHU_wm_list(l,1) = 3;
            elseif ismember(this_wm_roi,paralimbic_wm) 
                JHU_wm_list(l,1) = 4;
            elseif ismember(this_wm_roi,commissural)
                JHU_wm_list(l,1) = 5;
            else 
                JHU_wm_list(l,1) = 0;
            end
        end
%         wm_coords = elec_mni_coords(find(type_ind==-1),:);
%         wm_names = ch_names(find(type_ind==-1),:);
%         [mni_coords, wm_roi, NN_flag] = nifti_values(wm_coords,'JHU_MNI_SS_WMPM_Type-I_to_MNI_brain.nii');
%         
        all_wm_roi = zeros(length(type_ind),1);
        all_wm_roi(find(usable_inds),1) = JHU_wm_list;
        all_node_roi = zeros(length(type_ind),1);
        all_node_roi(find(usable_inds),1) = node_roi;
        all_lobe_list = zeros(length(type_ind),1);
        all_lobe_list(find(usable_inds),1) = lobe_list;
        
        % now we must make sure we only keep wm roi which are definitely wm
        % based on the segmentation
        all_non_wm = [find(type_ind==1);find(type_ind==0)];
        all_wm_roi(all_non_wm) = 0;
        
        patient_loc(polarity).session(session).coords = elec_mni_coords;
        patient_loc(polarity).session(session).names = final_ch_names;
        patient_loc(polarity).session(session).type = type_ind;
        patient_loc(polarity).session(session).roi = all_node_roi;
        patient_loc(polarity).session(session).lobes = all_lobe_list;
        patient_loc(polarity).session(session).t1 = this_T1_regions;
        patient_loc(polarity).session(session).seg = gm_seg_inds;
        patient_loc(polarity).session(session).eve = EVE_wm_only';
        patient_loc(polarity).session(session).seg_ID = gm_seg_label_ID;
        patient_loc(polarity).session(session).wm = all_wm_roi;
        patient_loc(polarity).session(session)
        
       
        % quadrant based analysis
        if polarity==1
            quadrant_ch_names = gchQ((ssIDq==session),:);
            this_session_quadrants = qexp_quad(ssIDq==session);
        else
            quadrant_ch_names = gchP((ssIDp==session),1);
            this_session_quadrants = pow_quad(ssIDp==session);
        end
        
        num_ch = length(quadrant_ch_names);
        for e = 1:num_ch
           try node_type(e) = type_ind(strcmp(patient_loc(polarity).session(session).names(:,1),quadrant_ch_names(e)));
           catch ME
               node_type(e) = 0;
           end
           
           try seg_type(e) = gm_seg_inds(strcmp(patient_loc(polarity).session(session).names(:,1),quadrant_ch_names(e)));
           catch ME
               seg_type(e) = 0;
           end
           
           try wm_type(e) = all_wm_roi(strcmp(patient_loc(polarity).session(session).names(:,1),quadrant_ch_names(e)));
           catch ME
               wm_type(e) = 0;
           end
        end
        
        if polarity==2
        quad_nodes_p = [quad_nodes_p;node_type'];
        quad_seg_p = [quad_seg_p;seg_type'];
        quad_wm_p = [quad_wm_p;wm_type'];
        else
        quad_nodes_q = [quad_nodes_q;node_type'];
        quad_seg_q = [quad_seg_q;seg_type'];
        quad_wm_q = [quad_wm_q;wm_type'];
        end
        
        all_T1_roi = [all_T1_roi;this_T1_regions];
    end

    
    
end
% region list
% Frontal WM
% Temporal WM

%% do quadrant analysis
plot_x_label = {'Front. ctx.','Temp-Par. ctx.','Paralimbic GM',...
                'Thal.-cort. WM.','Front.-assoc. WM','Temp.-assoc. WM',...
                'Paralimbic WM','Commissural WM','Background'};

code_regions_pow = quad_seg_p.*(quad_nodes_p==1)+(quad_wm_p+3).*(quad_nodes_p==-1);
code_regions_qexp = quad_seg_q.*(quad_nodes_q==1)+(quad_wm_q+3).*(quad_nodes_q==-1);

for i = [1,3]
    for j = [1:8]
        
        qexp_quad_plot(i,j) = sum((qexp_quad==i).*(code_regions_qexp==j));
        pow_quad_plot(i,j) = sum((pow_quad==i).*(code_regions_pow==j));
    end
end

qexp_quad_plot(2,:) = [];
pow_quad_plot(2,:) = [];

all_quad_data(1,:) = [qexp_quad_plot(1,:)-qexp_quad_plot(2,:)];
all_quad_data(2,:) = [pow_quad_plot(1,:)-pow_quad_plot(2,:)];


chi2Tests(qexp_quad_plot)
chi2Tests(pow_quad_plot)
figure(1);clf;
bar(all_quad_data)
legend(plot_x_label)
xticks([1:2])
xticklabels({'Communicability','Power'})
ylabel('Node count')
title('Anatomical distribution difference: quadrant 1 - quadrant 3')


%%

figure(1);clf;
subplot(1,2,1)
bar(qexp_quad_plot)
legend(plot_x_label)
xticks([1:2])
xticklabels({'Quadrant 1','Quadrant 3'})
ylabel('Node count')
ylim([0 80])
title(sprintf('Qexp, Chi-square p = %d',chi2Tests(qexp_quad_plot)))
subplot(1,2,2)
bar(pow_quad_plot)
legend(plot_x_label)
xticks([1:2])
xticklabels({'Quadrant 1','Quadrant 3'})
ylabel('Node count')
ylim([0 70])
title(sprintf('Power, Chi-square p = %f',chi2Tests(pow_quad_plot)))
figure(2);clf;
subplot(1,2,1)
bar(qexp_quad_plot)
legend('Frontal WM','Frontal GM','Temporal WM','Temporal GM')
xticks([1:4])
xticklabels({'Quadrant 1','Quadrant 2','Quadrant 3','Quadrant 4'})
ylabel('Node count')
title(sprintf('Qexp, Chi-square p = %d',chi2Tests(qexp_quad_plot)))
subplot(1,2,2)
bar(pow_quad_plot)
legend('Frontal WM','Frontal GM','Temporal WM','Temporal GM')
xticks([1:4])
xticklabels({'Quadrant 1','Quadrant 2','Quadrant 3','Quadrant 4'})
ylabel('Node count')
title(sprintf('Power, Chi-square p = %f',chi2Tests(pow_quad_plot)))

%% Do anatomical analysis

for r = 1:8
    which_subregions = all_sub_roi(r).data;
    all_region_counts(r).data = zeros(length(which_subregions),2,8);
end

plot_x_label = {'Front. ctx.','Temp-Par. ctx.','Paralimbic GM',...
                'Thal.-cort. WM.','Front.-assoc. WM','Temp.-assoc. WM',...
                'Paralimbic WM','Commissural WM','Background'};

qexp_vs_pow = {'qexp','pow'};
threshold = [0.05, 1.1, 1.1];
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];
color5 = [103 55 155]/255;
color6 = [78 172 91]/255;
load color_bar
folder_list = {'p_005','all','min_max'};
freq_bands = {'Alpha-theta','Beta','Low-gamma','High-gamma'};

z = 0;
yy = 0;
% do for qexp and pow
for metric = 1:2
    
    if metric==1
        load graphRT_sIall_500mspre
    else 
        load powRT_sIall_500mspre
    end
    
    clear EVE_all_freq
    
    % do for frequency bands 1:4
    for freq = 1:4
        z = z+1;
        
        
        % do for thresholds 0.05 and 0.1
        for thres = 1
             clear patient_roi
             clear mean_bval_lobe
             clear EVE_100_elec
             
             all_non_sig = [];
             pt_non_sig = [];
             seg_sig_coords = [];
             wm_sig_coords = [];
             seg_sig_all = [];
             wm_sig_all = [];
             background_sig_coords = [];
             background_pval_all = [];
             seg_pval_all = [];
             wm_pval_all = [];
             
             
             for r = 1:9
                seg_EVE_roi(r).data = [];
             end
             
             patient_roi = NaN*zeros(1,27);
             seg_bvals = NaN*zeros(27,3);
             seg_counts = NaN*zeros(27,9);
             EVE_100_elec = NaN*zeros(27,12);
             %clear mean_bval
             patient_roi_min = zeros(177,27);
             patient_roi_max = zeros(177,27);
            
            % do for sessions 1:27
            for session = 1:27
                %session

             
                % Get all b and p values
                if strcmp(qexp_vs_pow{metric},'qexp')
                    all_b_val = LTqexp{session}(:,freq,1);
                    all_p_val = LTqexp{session}(:,freq,4);
                else
                    all_b_val = LTpow{session}(:,freq,1);
                    all_p_val = LTpow{session}(:,freq,4);
                end

                % get significant channels based on p value threshold
                sig_ch = find(all_p_val<threshold(thres));
                non_sig_ch = find(all_p_val>=.05);
                ch_names = NFstruct(session).gchlbl;
                
                % *** still need to do channel parsing ***
                
                % check if there are any significant channels
                if isempty(sig_ch)
                    fprintf('skipping session b/c no significance\n')
                    pt_non_sig = [pt_non_sig;mean(all_b_val)];
                    all_non_sig = [all_non_sig;all_b_val];
                else
                    % extract significant channel p values, p values, and names
                    sig_p_val = all_p_val(sig_ch);
                    sig_b_val = all_b_val(sig_ch);
                    sig_ch_names = ch_names(sig_ch);
                    
                    non_sig_b_val = all_b_val(non_sig_ch);
                    pt_non_sig = [pt_non_sig;nanmean(non_sig_b_val)];
                    all_non_sig = [all_non_sig;non_sig_b_val];
                    
                    patient_coord_mapping = zeros(length(sig_ch_names),1);
                    for ch = 1:length(sig_ch_names)
                        % get this individual channel name
                        this_ch = sig_ch_names{ch,:};
                        % get coordinate index
                        coord_index = find(strcmp(patient_loc(metric).session(session).names(:,1),this_ch));
                        
                        try patient_coord_mapping(ch) = coord_index;
                        catch ME
                            fprintf('cannot find electrode %s\n', this_ch)
                            
                        end
                    end
                    
                    [ind_max] =  (sig_b_val==max(sig_b_val));
                    [ind_min] =  (sig_b_val==min(sig_b_val));
                    
                    % need csv files for all, <0.05, and min/max
                    
                    sig_roi = patient_loc(metric).session(session).roi(patient_coord_mapping(patient_coord_mapping~=0));
                    sig_type = patient_loc(metric).session(session).type(patient_coord_mapping(patient_coord_mapping~=0));
                    sig_lobe = patient_loc(metric).session(session).lobes(patient_coord_mapping(patient_coord_mapping~=0));
                    sig_seg = patient_loc(metric).session(session).seg(patient_coord_mapping(patient_coord_mapping~=0));
                    sig_wm = patient_loc(metric).session(session).wm(patient_coord_mapping(patient_coord_mapping~=0));
                    sig_mni = patient_loc(metric).session(session).coords(patient_coord_mapping(patient_coord_mapping~=0),:);
                    %sig_eve = patient_loc(metric).session(session).eve(patient_coord_mapping(patient_coord_mapping~=0),:);
                    sig_seg_ID = patient_loc(metric).session(session).seg_ID(patient_coord_mapping(patient_coord_mapping~=0),:);
                    proc_pval = sig_p_val(patient_coord_mapping~=0);
                    proc_bval = sig_b_val(patient_coord_mapping~=0);
                    
                    % do GM/WM analysis
                    mean_bval(session,(2*(freq-1)+1)) = mean(proc_bval(sig_type==-1));
                    mean_bval(session,(2*freq)) = mean(proc_bval(sig_type==1));
                    
                    if thres==1
                    frac_sig(session,(2*(freq-1)+1)) = sum(sig_type==-1)./sum(patient_loc(metric).session(session).type==-1);
                    frac_sig(session,(2*freq)) = sum(sig_type==1)./sum(patient_loc(metric).session(session).type==1);
                    end
                    
                    % GM/WM + lobe analysis
                    mean_bval_lobe(session,1) = mean(proc_bval([sig_type==-1]&[sig_lobe==1]));
                    mean_bval_lobe(session,2) = mean(proc_bval([sig_type==1]&[sig_lobe==1]));
                    mean_bval_lobe(session,3) = mean(proc_bval([sig_type==-1]&[sig_lobe==2]));
                    mean_bval_lobe(session,4) = mean(proc_bval([sig_type==1]&[sig_lobe==2]));
                    
                    d = 0;
                    for r = [0 4 12 19 20 27 37 68 69 70 82 83]
                        d = d+1;
                        if r==0
                            background_vals = sig_b_val(sig_roi==r);
                            EVE_100_elec(session,d) = nanmean(background_vals);
                        else
                            left_side = sig_b_val(sig_roi==r);
                            right_side = sig_b_val(sig_roi==(r+88));
                            EVE_100_elec(session,d) = nanmean([left_side;right_side]);
                        end
                    end
                    
                    % go through GM segmentation aggregate regions
                    for r = 1:3
                        seg_bvals(session,r) = nanmean(sig_b_val(sig_seg==r));
                        seg_counts(session,r) = length(sig_b_val(sig_seg==r));
                        seg_EVE_roi(r).data = [seg_EVE_roi(r).data;sig_b_val(sig_seg==r)];
                        seg_sig_all = [seg_sig_all; sig_b_val(sig_seg==r)];
                        seg_pval_all = [seg_pval_all; sig_p_val(sig_seg==r)];
                        seg_sig_coords = [seg_sig_coords; sig_mni((sig_seg==r),:)];

                    end
                    
                    % go through WM EVE aggregate regions
                    for r = 1:5
                        EVE_WM(session,r) = nanmean(sig_b_val(sig_wm==r));
                        seg_counts(session,(r+3)) = length(sig_b_val(sig_wm==r));
                        seg_EVE_roi(r+3).data = [seg_EVE_roi(r+3).data;sig_b_val(sig_wm==r)];
                        wm_sig_all = [wm_sig_all; sig_b_val(sig_wm==r)];
                        wm_pval_all = [wm_pval_all; sig_p_val(sig_wm==r)];
                        wm_sig_coords = [wm_sig_coords; sig_mni((sig_wm==r),:)];
                    end
                    
%                     for r = 1:8
%                         which_subregions = all_sub_roi(r).data;
%                         for qq = 1:length(which_subregions)
%                             if r<4 
%                                 sub_seg_high = find((sig_seg_ID==which_subregions(qq)) .*(proc_bval >0));
%                                 sub_seg_low = find((sig_seg_ID==which_subregions(qq)) .*(proc_bval <0));
%                             else
%                                 sub_seg_high = find((sig_eve==which_subregions(qq)) .*(proc_bval >0));
%                                 sub_seg_low = find((sig_eve==which_subregions(qq)) .*(proc_bval <0));
%                             end
%                             
%                             all_region_counts(r).data(qq,1,z) = all_region_counts(r).data(qq,1,z)+length(sub_seg_high)
%                             all_region_counts(r).data(qq,2,z) = all_region_counts(r).data(qq,2,z)+length(sub_seg_low)
%                         end
%                         
%                     end
                    
                    
                    seg_counts(session,9) = length(sig_b_val(sig_type==0));
                    
                    pt_background(session,1) = nanmean(sig_b_val(sig_type==0));
                    seg_EVE_roi(9).data = [seg_EVE_roi(9).data;sig_b_val(sig_type==0)];
                    background_sig_coords = [background_sig_coords; sig_mni((sig_type==0),:)];
                    background_pval_all = [background_pval_all;sig_p_val(sig_type==0)];
                    
                    % deterine how many nodes are in each roi
                    if thres<3
                        for r = 0:176
                            row_ind = r+1;
                            patient_roi(row_ind,session) = sum(sig_roi==r);
                        end
                        
                    else
                        % get the indices (roi) for the minimum value node
                        row_ind_min = patient_loc(metric).session(session).roi(ind_min)+1;
                        patient_roi_min(row_ind_min,session) = 1;
                        
                        % get the indices (roi) for the max value node
                        row_ind_max = patient_loc(metric).session(session).roi(ind_max)+1;
                        patient_roi_max(row_ind_max,session) = 1;
                        
                    end
                    
                end
                
                
            end
            
%             figure(z);clf;
%             hold on
%             imagesc(seg_counts)
%             xticks([1:9])
%             xticklabels(plot_x_label)
%             yticks([1:27])
%             yticklabels(session_list)
%             title(sprintf('Metric = %s Frequency = %s node counts',qexp_vs_pow{metric},freq_bands{freq}))
%             colorbar
%             set(gcf,'position',[100,100,1200,900])
%             %saveas(gcf,sprintf('output/Final/dec28_plots/region_counts/metric_%s_frequency_%s_seg_EVE.jpg',qexp_vs_pow{metric},freq_bands{freq}))
%             hold off
            
            
            feat_nums = sum(seg_counts(:,1:9),2);
            feat_nums(isnan(feat_nums)) = 0;
            feat_gen(1:27,z) = feat_nums;
            
            
            all_non_sig_struct(z).data = all_non_sig;
            pt_non_sig_struct(z).data = pt_non_sig;
            
            EVE_all_freq(:,:,freq) = EVE_100_elec;
            
            data_table = [seg_bvals,EVE_WM,pt_background];
            
            for r = 1:8
                all_median_table(z,r) = nanmedian(seg_EVE_roi(r).data);  
                all_pval(z,r) = ranksum(seg_EVE_roi(r).data,all_non_sig_struct(z).data);
            end
            
            figure(z);clf
            hold on
            qqq= 0;
            all_or_pt=0;
            if all_or_pt==0
            a1 = scatter(ones(1,length(all_non_sig_struct(z).data)),all_non_sig_struct(z).data,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
            plot([0.75 1.25],[nanmedian(all_non_sig_struct(z).data) nanmedian(all_non_sig_struct(z).data)],'k-','LineWidth',2)
            for r = 1:9
                if r<4
                    which_color = color2;
                elseif r<9
                    which_color = color5;
                else
                    which_color = color1;
                end
                try p99(z).data(r) = ranksum(seg_EVE_roi(r).data,all_non_sig_struct(z).data);
                        qqq = r+2;
                        scatter((r+2)*ones(1,length(seg_EVE_roi(r).data)),seg_EVE_roi(r).data,'MarkerEdgeColor',which_color,'MarkerFaceColor',which_color,'jitter','on')
                        plot([r+1.75 r+2.25],[nanmedian(seg_EVE_roi(r).data) nanmedian(seg_EVE_roi(r).data)],'k-','LineWidth',2)

                catch ME
                    p99(z).data(r) = 1;
                end
            end
            for r = 1:11
            try max_regions(r,1) = nanmax(seg_EVE_roi(r).data);
            catch ME
                max_regions(r,1) = 0;
            end
            end
            y_val = nanmax([max_regions;all_non_sig_struct(z).data(:)]);
            num_sig = sum(p99(z).data<(0.05/64));
            which_sig = find(p99(z).data<(0.05/64));
            for j = 1:num_sig
                stop_sig = which_sig(j)+2;
            plot([1,stop_sig], [y_val*1.1.^j,y_val*1.1.^j], '-k', 'LineWidth',2)
            if p99(z).data(stop_sig-2)<(0.001./64)
            plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), '+k')
            elseif p99(z).data(stop_sig-2)<(0.01./64)
                plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), 'dk')
            elseif p99(z).data(stop_sig-2)<(0.05./64)
                plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), '*k')
%             elseif p99(z).data(stop_sig-2)<0.05
%                 plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), 'xk')
            end
            end
            hold off
%             try ylim([y_val*(1.1.^j+0.1)])
%             catch ME
%             end
%             legend([a5,a6],'White matter','Gray matter','Location','SouthEast')
%             xticks([1,3.5,5.5])
%             ylabel('Regression slope')
             xticks([1,3:qqq])
             xticklabels(['Non-significant',plot_x_label])
            set(gcf,'position',[100,100,1333,500])
            title(sprintf('Metric = %s Frequency = %s lobar ROI',qexp_vs_pow{metric},freq_bands{freq}))
            %saveas(gcf,sprintf('output/10jan2021/inst0/metric_%s_frequency_%s_seg_EVE.svg',qexp_vs_pow{metric},freq_bands{freq}))
            
            
            elseif all_or_pt==1
                % pt level
            a1 = scatter(ones(1,length(pt_non_sig_struct(z).data)),pt_non_sig_struct(z).data,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
            plot([0.75 1.25],[nanmedian(pt_non_sig_struct(z).data) nanmedian(pt_non_sig_struct(z).data)],'k-','LineWidth',2)
            for r = 1:size(data_table,2)
                try p99(z).data(r) = ranksum(data_table(:,r),all_non_sig_struct(z).data);
                    if r<4
                    which_color = color2;
                    elseif r<9
                    which_color = color5;
                    else
                    which_color = color1;
                    end
                    if sum(~isnan(data_table(:,r)))>=2
                        qqq = r+2;
                        scatter((r+2)*ones(1,27),data_table(:,r),'MarkerEdgeColor',which_color,'MarkerFaceColor',which_color,'jitter','on')
                        plot([r+1.75 r+2.25],[nanmedian(data_table(:,r)) nanmedian(data_table(:,r))],'k-','LineWidth',2)
                    else 
                        p99(z).data(r) = 1;
                        
                    end
                catch ME
                    p99(z).data(r) = 1;
                end
            end
            
            y_val = nanmax([data_table(:);all_non_sig_struct(z).data(:)]);
            num_sig = sum(p99(z).data<0.05);
            which_sig = find(p99(z).data<0.05);
            for j = 1:num_sig
                stop_sig = which_sig(j)+2;
            plot([1,stop_sig], [y_val*1.1.^j,y_val*1.1.^j], '-k', 'LineWidth',2)
            if p99(z).data(stop_sig-2)<(0.001./64)
            plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), '+k')
            elseif p99(z).data(stop_sig-2)<(0.01./64)
                plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), 'dk')
            elseif p99(z).data(stop_sig-2)<(0.05./64)
                plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), '*k')
            elseif p99(z).data(stop_sig-2)<0.05
                plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), 'xk')
            end
            end
            hold off
%             try ylim([y_val*(1.1.^j+0.1)])
%             catch ME
%             end
%             legend([a5,a6],'White matter','Gray matter','Location','SouthEast')
%             xticks([1,3.5,5.5])
%             ylabel('Regression slope')
             xticks([1,3:qqq])
            xticklabels(['Non-significant',plot_x_label])
            set(gcf,'position',[100,100,1333,500])
            title(sprintf('Metric = %s Frequency = %s lobar ROI',qexp_vs_pow{metric},freq_bands{freq}))
            %saveas(gcf,sprintf('output/Final/dec21_plots/seg_GM_EVE_pt/metric_%s_frequency_%s_seg_EVE.jpg',qexp_vs_pow{metric},freq_bands{freq}))
            
                
            end
             
             
             for r = 1:4
                p100(z).data(r) = ranksum(mean_bval_lobe(:,r),all_non_sig_struct(z).data);
            end
            

            if thres <3
                % create table
                patient_roi_bilat = zeros(89,27);
                patient_roi_bilat(1,:) = patient_roi(1,:);
                patient_roi_bilat([2:89],:) =  patient_roi_bilat([2:89],:) + patient_roi(2:89,:)+patient_roi(90:177,:);
                region_table = create_region_table(EVE_roi_list_bilat,patient_roi_bilat);
                % write table
                which_folder = folder_list{thres};
%                 writetable(region_table,sprintf('output/Final/CSV_files/%s/%s_counts_freq_%d_%f_threshold.csv',...
%                     which_folder,qexp_vs_pow{metric},freq,threshold(thres))); 
            else
                 % create the tables
                patient_roi_bilat_min = zeros(89,27);
                patient_roi_bilat_min(1,:) = patient_roi_min(1,:);
                patient_roi_bilat_min([2:89],:) =  patient_roi_bilat_min([2:89],:) + patient_roi_min(2:89,:)+patient_roi_min(90:177,:);
                
                patient_roi_bilat_max = zeros(89,27);
                patient_roi_bilat_max(1,:) = patient_roi_max(1,:);
                patient_roi_bilat_max([2:89],:) =  patient_roi_bilat_max([2:89],:) + patient_roi_max(2:89,:)+patient_roi_max(90:177,:);

                region_table_min = create_region_table(EVE_roi_list_bilat,patient_roi_bilat_min);
                region_table_max = create_region_table(EVE_roi_list_bilat,patient_roi_bilat_max);

                which_folder = folder_list{thres};

                % write the tables
%                 writetable(region_table_min,sprintf('output/Final/CSV_files/%s/%s_min_counts_freq_%d_%f_threshold.csv',...
%                     which_folder,qexp_vs_pow{metric},freq,threshold(thres))); 
% 
%                 writetable(region_table_max,sprintf('output/Final/CSV_files/%s/%s_max_counts_freq_%d_%f_threshold.csv',...
%                     which_folder,qexp_vs_pow{metric},freq,threshold(thres))); 
                    
            end

        end
        
        % do rendering of three categories:
        % (1) all nodes in three groups
%         
        size(seg_sig_all)
        size(seg_sig_coords)
        
        size(wm_sig_all)
        size(wm_sig_coords)
        
        size(seg_EVE_roi(9).data)
        size(background_sig_coords)
        
        % gm
        final_elec_matrix = [seg_sig_coords,sign(-1.*seg_sig_all),1./seg_pval_all];
        %dlmwrite('output/1Jan2021/nodal_plots/render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
        %BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/1Jan2021/nodal_plots/render_elecs.node','dec28_option.mat',sprintf('output/1Jan2021/nodal_plots/GM_nodal_%d_freq_%d.jpg',metric,freq))
         
        % wm
        final_elec_matrix = [wm_sig_coords,sign(-1.*wm_sig_all),1./wm_pval_all];
        %dlmwrite('output/1Jan2021/nodal_plots/render_elecs.node',final_elec_matrix,'delimiter',' ','precision',5)
        %BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','output/1Jan2021/nodal_plots/render_elecs.node','dec28_option.mat',sprintf('output/1Jan2021/nodal_plots/WM_nodal_%d_freq_%d.jpg',metric,freq))

        %delete output/1Jan2021/nodal_plots/render_elecs.node
        
        % do rendering
        % mean value of each ROI
        for r = 1:8
            region_values(r) = nanmean(seg_EVE_roi(r).data);
            region_pvals(r) = p99(z).data(r);
        end

       min_val = min(region_values);
       max_val = max(region_values);

       if abs(min_val)>abs(max_val)
           entry = floor(region_values./abs(min_val)*50)+51;
       else
           entry = floor(region_values./abs(max_val)*50)+50;
       end

       if sum(isnan(entry))
           entry(isnan(entry)) =50;
       end
      
       for r = 1:8
           if region_pvals(r)>=(0.05./64)
               entry(r) = 50;
           end
       end
           
       color_map_regions = color_bar(entry,:);
           
       cmap_GM_raw = color_map_regions(1:3,:);
       cmap_WM_raw = color_map_regions(4:8,:);
       
       % find number of wm 1
       l1 = length(thalamo_cortical);
       l2 = length(frontal_association);
       l3 = length(temporal_association);
       l4 = length(paralimbic_wm);
       l5 = length(commissural);
       
       cmap_WM = zeros((l1+l2+l3+l4+l5),3);
       for l = 1:l1
       cmap_WM(l,:)=cmap_WM_raw(1,:);
       end
       for l = (l1+1):(l1+l2)
       cmap_WM(l,:)=cmap_WM_raw(2,:);
       end
       for l = (l1+l2+1):(l1+l2+l3)
       cmap_WM(l,:)=cmap_WM_raw(3,:);
       end
       for l = (l1+l2+l3+1):(l1+l2+l3+l4)
       cmap_WM(l,:)=cmap_WM_raw(4,:);
       end
       for l = (l1+l2+l3+l4+1):(l1+l2+l3+l4+l5)
       cmap_WM(l,:)=cmap_WM_raw(5,:);
       end
       
       l6 = 12;
       l7 = 16;
       l8 = 14;
       
       cmap_GM = zeros(40,3);
       for l = 1:l6
       cmap_GM(l,:) = cmap_GM_raw(1,:);
       end
       for l = (l6+1):(l6+l7)
       cmap_GM(l,:) = cmap_GM_raw(2,:);
       end
       for l = (l6+l7+1):(l6+l7+l8)
       cmap_GM(l,:) = cmap_GM_raw(3,:);
       end


        
           %now we must map colors into composite regions
           dlmwrite(sprintf('output/1Jan2021/surface_plots/cmap_GM_metric_%f_freq_%f.txt',metric,freq),cmap_GM,'delimiter',' ','precision',5)
           dlmwrite(sprintf('output/1Jan2021/surface_plots/cmap_WM_metric_%f_freq_%f.txt',metric,freq),cmap_WM,'delimiter',' ','precision',5)
           
           
           %BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','JHU48.nii','JHU_surface_rendering.mat',sprintf('WM_surface_metric_%d_freq_%d.jpg',metric,freq))
           %BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii','Yeo_surface_rendering.mat',sprintf('GM_surface_metric_%d_freq_%d.jpg',metric,freq))
         
        
    end
    
    

%     
%     %boxplot([mean_bval],'Labels',{'WM','GM','WM','GM','WM','GM','WM','GM'})
%     c_pos = multcompare(stats,[],'off');
%     start_sig = c_pos(find(c_pos(:,6)<0.05),1);
%     stop_sig = c_pos(find(c_pos(:,6)<0.05),2);
% 
%     ylabel('Mean patient regression slope')
%     title(sprintf('%s metric, Kruskal-Wallis test p = %f',qexp_vs_pow{metric},p));
%     if ~isempty(start_sig)
%         num_sig = length(start_sig);
%         for j = 1:num_sig
%         y_val = nanmax([mean_bval(:,start_sig(j));mean_bval(:,stop_sig(j))]);
%         plot([start_sig(j),stop_sig(j)], [y_val*1.1.^j,y_val*1.1.^j], '-k', 'LineWidth',2)
%         plot(mean([start_sig(j),stop_sig(j)]), y_val*(1.1.^j+0.05), '*k')
%         end
%     end
%     ylim([-600 1100])
%     xticks([1:8])
%     xticklabels({'WM','GM','WM','GM','WM','GM','WM','GM'})
%     legend([a1,a3,a5,a7],'Alpha/Theta','Beta','Low-gamma','High-gamma','Location','SouthWest')
%     hold off
            
end

non_sigs = find(all_pval>(0.05));

all_median_table(non_sigs) = NaN;

figure(9);clf;
subplot(1,2,1)
hold on
imagesc(all_median_table(1:4,:)','AlphaData',~isnan(all_median_table(1:4,:)'))
colormap(color_bar)
caxis([-200 200])
ylim([0.5 8.5])
xlim([0.5 4.5])
colorbar
hold off
subplot(1,2,2)
hold on
imagesc(all_median_table(5:8,:)','AlphaData',~isnan(all_median_table(5:8,:)'))
colormap(color_bar)
caxis([-350 350])
ylim([0.5 8.5])
xlim([0.5 4.5])
colorbar
hold off
%%
        for r = 1:8
            for z = 1:8
                this_roi_data = all_region_counts(r).data(:,:,z);
                
                this_roi_data(find(sum(this_roi_data,2)==0),:) = [];
                try subregion_pvals(r,z) = chi2Tests(this_roi_data)
                catch ME
                    subregion_pvals(r,z) = NaN;
                end
            end
        end

% save('all_non_sig_struct.mat','all_non_sig_struct')
% save('pt_non_sig_struct.mat','pt_non_sig_struct')

