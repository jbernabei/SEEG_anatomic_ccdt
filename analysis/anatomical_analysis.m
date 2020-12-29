%% anatomical_analysis.m
% John Bernabei
% June 2020
% Center for Neuroengineering & Therapeutics

%% set up workspace
clear all

% load required data
load nodeIDsForJohn

load all_regions_raw

% Region list
temporal_list = [1:16,30:31,82:83];
frontal_list = [20:21, 24:29, 50:59, 68:73, 76:81];
occipital_list = [64:67, 22:23];
parietal_list = [32:33, 60:63];

all_WM_data =readtable('/Users/jbernabei/Downloads/Atlas_labels/JHU48.txt');
%all_GM_data =readtable('/Users/jbernabei/Downloads/Yeo_JNeurophysiol11_MNI152/Yeo2011_7Networks_ColorLUT.txt');

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
%% Data processing figure for all patients

all_vs_minmax = 'all'; % 'all' if we want all nodes, other if we want min & max only
qexp_vs_pow = 'qexp'; % 'qexp' if we want qexp, 'pow' if we want power
threshold = 'high'; % 'low' if we wand <0.05, 'high' if we want <0.10

% load separate label files depending on metric
if strcmp(qexp_vs_pow,'qexp')
    load graphChLbls
else 
    load powChLbls
end

mean_bval_temp_GM_all_freq = [];
mean_bval_temp_WM_all_freq = [];
mean_bval_front_GM_all_freq = [];
mean_bval_front_WM_all_freq = [];

observed_frac_front_GM_pos_all_freq = [];
observed_frac_front_WM_pos_all_freq = [];
observed_frac_front_GM_neg_all_freq = [];
observed_frac_front_WM_neg_all_freq = [];

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