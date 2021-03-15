pfnm = 'powRT_sIall_500mspre';
nfnm = 'plvRT_sIall_100mspre';
gfnm = 'graphRT_sIall_500mspre';

ploton=0;
nodeType = 0; %0 = all feature nodes, 1 = min/max slope, 2 = min P

load(['/Users/jbernabei/Documents/PhD_Research/Neurosurgery_SEEG/SEEG_anatomic_ccdt/' pfnm '.mat'],'NFstruct', 'LTpow');
NFstructP = NFstruct;
load(['/Users/jbernabei/Documents/PhD_Research/Neurosurgery_SEEG/SEEG_anatomic_ccdt/' nfnm '.mat']);
NFstructPLV = NFstruct;
LTplv = LTqexp;
load(['/Users/jbernabei/Documents/PhD_Research/Neurosurgery_SEEG/SEEG_anatomic_ccdt/' gfnm '.mat']);

tS=[];
tSP = [];
fQv = [];
sQv = [];
fPv = [];
sPv = [];
fQvPLV = [];
sQvPLV = [];
fQvns = [];
sQvns = [];
fPvns = [];
sPvns = [];
gchQ = [];
gchP = [];
ssIDq = [];
ssIDp = [];
fIDp = [];
fIDq = [];
ssIDplv = [];
fIDplv = [];
nodeIDp1 = [];
nodeIDp2 = [];
nodeIDq = [];
mniQroi = [];
mniProi = [];
yeoQroi = [];
yeoProi = [];
db = CCDTdatabase;
agmniQ = [];
agmniP = [];
agyeoQ = [];
agyeoP = [];
y=0;
z=0;
w=0;
u=0;
uu=0;
%%
for x=1:27
    cgQx = NFstruct(x).gchlbl;
    cgPx = NFstructP(x).gchlbl;
    cgPLVx = NFstructPLV(x).gchlbl;
    fid = fopen(['/Users/jbernabei/Documents/PhD_Research/Neurosurgery_SEEG/improved_localization/new_' db{x,1} '_localization.csv']);
    J = textscan(fid,'%s %s %f %f %f %s %s', 'delimiter',',');
    fclose(fid);
    cch = J{1};
    t1 = J{2};
    mni = J{6};
    yeo = J{7};
    uchQ = ismember(cch,cgQx);
    uchP = ismember(cch,cgPx(:,1));
    uchPLV = ismember(cch,cgPLVx);

    cmniQ = mni(uchQ);
    cmniP = mni(uchP);
    cyeoQ = yeo(uchQ);
    cyeoP = yeo(uchP);
    
    if isequal(length(cmniQ),length(cgQx))
        z=z+1;
        goodSubj{z,1} = db{x,1};
        goodSubj{z,2} = cgQx(~ismember(cgQx,cch));
    else
        y=y+1;
        badSubj{y,1} = db{x,1};
        badSubj{y,2} = cgQx(~ismember(cgQx,cch));
    end
    
    for j=1:4
        czQ = zscore(NFstruct(x).fbands(j).NFqexp_nodal,1);
        czP = zscore(NFstructP(x).fbands(j).NFpow,1);
        czPLV = zscore(NFstructPLV(x).fbands(j).NFqexp_nodal,1);

        if ~nodeType %all feature selected nodes
            cLT = LTqexp{x}(:,j,4);
            cLTP = LTpow{x}(:,j,4);
            cLTplv = LTplv{x}(:,j,4);

            if sum(cLT<=0.05)>0
                cNodesQ = find(cLT<=0.05);
            else
                [p,cb] = min(cLT);
                u=u+1;
                minTrack(1) = u;
                pTrackQ(u) = p;
                cNodesQ = cb;
            end
            if sum(cLTP<=0.05)>0
                cNodesP = find(cLTP<=0.05);
            else
                [p,cb] = min(cLTP);
                w=w+1;
                minTrack(2) = w;
                pTrackP(w) = p;
                cNodesP = cb;
            end
            if sum(cLTplv<=0.05)>0
                cNodesPLV = find(cLTplv<=0.05);
            else
                [p,cb] = min(cLTplv);
                uu=uu+1;
                minTrack(3) = uu;
                pTrackPLV(uu) = p;
                cNodesPLV = cb;
            end
        elseif nodeType ==1 %min/max slope
            cLT = LTqexp{x}(:,j,1);
            cLTP = LTpow{x}(:,j,1);
            cLTplv = LTplv{x}(:,j,1);
            [~,cb] = min(cLT);
            [~, cb1] = max(cLT);
            cNodesQ = [cb cb1];
            [~,cb] = min(cLTP);
            [~, cb1] = max(cLTP);
            cNodesP = [cb cb1];
            [~,cb] = min(cLTplv);
            [~, cb1] = max(cLTplv);
            cNodesPLV = [cb cb1];
        else %single min p
            cLT = LTqexp{x}(:,j,4);
            cLTP = LTpow{x}(:,j,4);
            cLTplv = LTplv{x}(:,j,4);
            [~,cb] = min(cLT);
            cNodesQ = cb;
            [~,cb] = min(cLTP);
            cNodesP = cb;
            [~,cb] = min(cLTplv);
            cNodesPLV = cb;
        end
        cF = czQ(cNodesQ,behavStruct(x).ifast);
        cS = czQ(cNodesQ,behavStruct(x).islow);
        cgQ = cgQx(cNodesQ);
        
        cmniQroi = cmniQ(ismember(cch(uchQ),cgQ));
        cyeoQroi = cyeoQ(ismember(cch(uchQ),cgQ));
        cFns = czQ(~cNodesQ,behavStruct(x).ifast);
        cSns = czQ(~cNodesQ,behavStruct(x).islow);
        
        
        cFP = czP(cNodesP,behavStruct(x).ifast);
        cSP = czP(cNodesP,behavStruct(x).islow);
        cgP = cgPx(cNodesP,:);
        
        
        cFPns = czP(~cNodesP,behavStruct(x).ifast);
        cSPns = czP(~cNodesP,behavStruct(x).islow);
        cmniProi = cmniP(ismember(cch(uchP),cgP(:,1)));
        cyeoProi = cyeoP(ismember(cch(uchP),cgP(:,1)));
        
        cFplv = czPLV(cNodesPLV,behavStruct(x).ifast);
        cSplv = czPLV(cNodesPLV,behavStruct(x).islow);
        cgPLV = cgPLVx(cNodesPLV);
        
        sigNodeJkshtAnat{x} = cgQ(~ismember(cgQ,cch));
        
        
        [~,~,~,cstats] = ttest(cF,cS,'dim',2);
        [~,~,~,cstatsP] = ttest(cFP,cSP,'dim',2);
        cF = mean(cF,2);
        cS = mean(cS,2);
        cFP = mean(cFP,2);
        cSP = mean(cSP,2);
        cFplv = mean(cFplv,2);
        cSplv = mean(cSplv,2);
        cFns = mean(cFns,2);
        cSns = mean(cSns,2);
        cFPns = mean(cFPns,2);
        cSPns = mean(cSPns,2);
        cgP1 = cgP(:,1);
        cgP2 = cgP(:,2);
        csIDq = ones(length(cF),1).*x;
        csIDp = ones(length(cFP),1).*x;
        csIDplv = ones(length(cFplv),1).*x;
        cfq = ones(length(cF),1).*j;
        cfp = ones(length(cFP),1).*j;
        cfplv = ones(length(cFplv),1).*j;

        for i = 1:size(cgP,1)
            cnIDp1{i} = strcat([sprintf('%02d',x), num2str(j), cgP1{i}]);
            cnIDp2{i} = strcat([sprintf('%02d',x), num2str(j), cgP2{i}]);
        end
        for i = 1:size(cgQ,1)
            cnIDq{i} = strcat([sprintf('%02d',x), num2str(j), cgQ{i}]);
        end
        
        tS = [tS; cstats.tstat];
        tSP = [tSP; cstatsP.tstat];
        fQv = [fQv; cF];
        sQv = [sQv; cS];
        fPv = [fPv; cFP];
        sPv = [sPv; cSP];
        fQvPLV = [fQvPLV; cFplv];
        sQvPLV = [sQvPLV; cSplv];
        fQvns = [fQvns; cFns];
        sQvns = [sQvns; cSns];
        fPvns = [fPvns; cFPns];
        sPvns = [sPvns; cSPns];
        gchQ = [gchQ; cgQ];
        gchP = [gchP; cgP];
        agmniQ = [agmniQ; cmniQ];
        agmniP = [agmniP; cmniP];
        agyeoQ = [agyeoQ; cyeoQ];
        agyeoP = [agyeoP; cyeoP];
        mniQroi = [mniQroi; cmniQroi];
        mniProi = [mniProi; cmniProi];
        yeoQroi = [yeoQroi; cyeoQroi];
        yeoProi = [yeoProi; cyeoProi];
        ssIDq = [ssIDq; csIDq];
        ssIDp = [ssIDp; csIDp];
        ssIDplv = [ssIDplv; csIDplv];
        fIDp = [fIDp; cfp];
        fIDq = [fIDq; cfq];
        fIDplv = [fIDplv; cfplv];
        nodeIDp1 = [nodeIDp1; cnIDp1'];
        nodeIDp2 = [nodeIDp2; cnIDp2'];
        nodeIDq = [nodeIDq; cnIDq'];
        clear cnIDp* cnIDq
    end
end

%%
highNodesP = (fPv>0&sPv>0);
highNodesQ = (fQv>0&sQv>0);

hnNodeIDp1 = nodeIDp1(highNodesP);
hnNodeIDp2 = nodeIDp2(highNodesP);
hnNodeIDq = nodeIDq(highNodesQ);

c1 = ismember(nodeIDq,hnNodeIDp1);
c2 = ismember(nodeIDq,hnNodeIDp2);
c3 = c1+c2;
nodeShare_q_hp = logical(c3);
c1 = ismember(nodeIDq,nodeIDp1);
c2 = ismember(nodeIDq,nodeIDp2);
c3 = c1+c2;
nodeShare_q_p = logical(c3);

c1 = ismember(nodeIDp1,hnNodeIDq);
c2 = ismember(nodeIDp2,hnNodeIDq);
c3 = c1+c2;
nodeShare_p_hq = logical(c3);

c1 = ismember(nodeIDp1,nodeIDq);
c2 = ismember(nodeIDp2,nodeIDq);
c3 = c1+c2;
nodeShare_p_q = logical(c3);


clear c* i j x

% Separation
figure;
scatter(fQv,sQv,fIDq.*10,ssIDq)
title(['Feature Node Normalized Network Rank in Fast vs. Slow Trials: Qexp, nodeType ' num2str(nodeType)])
xlim([-2 2])
ylim([-2 2])
figure;
scatter(fPv,sPv,fIDp.*10,ssIDp)
title(['Feature Node Normalized Network Rank in Fast vs. Slow Trials: POW, nodeType ' num2str(nodeType)])
xlim([-2 2])
ylim([-2 2])
figure;
scatter(fQvPLV,sQvPLV,fIDplv.*10,ssIDplv)
title(['Feature Node Normalized Network Rank in Fast vs. Slow Trials: PLV, nodeType ' num2str(nodeType)])
xlim([-2 2])
ylim([-2 2])

cfQv = fQv(fIDq==1);
csQv = sQv(fIDq==1);
diagQ = abs(cfQv-csQv)/sqrt(2); %dist from diag
dQ = sqrt(cfQv.^2 + csQv.^2); %euclid dist from 0,0
rNpq = (abs(cfQv-csQv)/sqrt(2)) + (sqrt(cfQv.^2 + csQv.^2)); % distance from diagnoal + distance from 0,0
cfQv = fQv(fIDq==2);
csQv = sQv(fIDq==2);
diagQ2 = abs(cfQv-csQv)/sqrt(2); %dist from diag
dQ2 = sqrt(cfQv.^2 + csQv.^2); %euclid dist from 0,0
rNpq2 = (abs(cfQv-csQv)/sqrt(2)) + (sqrt(cfQv.^2 + csQv.^2)); % distance from diagnoal + distance from 0,0
cfQv = fQv(fIDq==3);
csQv = sQv(fIDq==3);
diagQ3 = abs(cfQv-csQv)/sqrt(2); %dist from diag
dQ3 = sqrt(cfQv.^2 + csQv.^2); %euclid dist from 0,0
rNpq3 = (abs(cfQv-csQv)/sqrt(2)) + (sqrt(cfQv.^2 + csQv.^2)); % di
cfQv = fQv(fIDq==4);
csQv = sQv(fIDq==4);
diagQ4 = abs(cfQv-csQv)/sqrt(2); %dist from diag
dQ4 = sqrt(cfQv.^2 + csQv.^2); %euclid dist from 0,0
rNpq4 = (abs(cfQv-csQv)/sqrt(2)) + (sqrt(cfQv.^2 + csQv.^2)); %


dQall = sqrt(fQv.^2 + sQv.^2);

dPall = sqrt(fPv.^2 + sPv.^2);
dPLVall = sqrt(fQvPLV.^2 + sQvPLV.^2);

diagQall = abs(fQv - sQv)./sqrt(2);
diagPLVall = abs(fQvPLV - sQvPLV)./sqrt(2);
diagPall = abs(fPv - sPv)./sqrt(2);


% figure;
% errorbar(1,mean(dQall),std(dQall))
% hold on
% errorbar(2,mean(dPLVall),std(dPLVall))
% errorbar(3,mean(dPall),std(dPall))
% bar([mean(dQall) mean(dPLVall) mean(dPall)])
% xticks(1:3)
% xticklabels({'Q' 'PLV' 'POW'})
% ylabel('Separation (Euclidean distance from 0,0)')
% [~,p,~,stats] = ttest2(dQall,dPLVall)
% [~,p,~,stats] = ttest2(dQall,dPall)

sem = std(diagQall)/sqrt(length(diagQall));
CI95 = tinv([0.025 0.975], length(diagQall)-1);                    
qCI95 = bsxfun(@times, sem, CI95(:)); 

sem = std(diagPLVall)/sqrt(length(diagPLVall));
CI95 = tinv([0.025 0.975], length(diagPLVall)-1);                    
plvCI95 = bsxfun(@times, sem, CI95(:)); 

sem = std(diagPall)/sqrt(length(diagPall));
CI95 = tinv([0.025 0.975], length(diagPall)-1);                    
pCI95 = bsxfun(@times, sem, CI95(:)); 

% figure;
% errorbar(1,mean(diagQall),qCI95(1),qCI95(2),'color',[0 0 0])
% hold on
% errorbar(2,mean(diagPLVall),plvCI95(1),plvCI95(2),'color',[0 0 0])
% errorbar(3,mean(diagPall),pCI95(1),pCI95(2),'color',[0 0 0])
% bar(1, mean(diagQall),'facecolor',[0.2 0.2 0.2])
% bar(2, mean(diagPLVall),'facecolor',[0.5 0.5 0.5])
% bar(3, mean(diagPall),'facecolor',[0.7 0.7 0.7])
% xticks(1:3)
% xticklabels({'Qexp' 'PLV' 'POW'})
% ylabel('Separation (Euclidean distance from diagonal +/- 95%CI)')
% [~,p,~,stats] = ttest2(diagQall,diagPLVall)
% [~,p,~,stats] = ttest2(diagQall,diagPall)
% title('Rank Separation of Feature Selected Nodes from Behavioral Equivalent Network Position')
% box off
% set(gca,'FontSize',14)
% set(gcf,'color','w')

figure;
errorbar(1,mean(diagQall),qCI95(1),qCI95(2),'color',[0 0 0])
hold on
errorbar(2,mean(diagPall),pCI95(1),pCI95(2),'color',[0 0 0])
bar(1, mean(diagQall),'facecolor',[0.2 0.2 0.2])
bar(2, mean(diagPall),'facecolor',[0.7 0.7 0.7])
xticks(1:2)
xticklabels({'Qexp' 'POW'})
ylabel('Separation (Euclidean distance from diagonal +/- 95%CI)')
[~,p,~,stats] = ttest2(diagQall,diagPall)
title('Feature Selected Nodes Distance from Behavioral Equivalent Diagonal Network Position')
box off
set(gca,'FontSize',14)
set(gcf,'color','w')

sem = std(dQ)/sqrt(length(dQ));
CI95 = tinv([0.025 0.975], length(dQ)-1);                    
dQCI95 = bsxfun(@times, sem, CI95(:)); 


sem = std(dQ2)/sqrt(length(dQ2));
CI95 = tinv([0.025 0.975], length(dQ2)-1);                    
dQ2CI95 = bsxfun(@times, sem, CI95(:)); 


sem = std(dQ3)/sqrt(length(dQ3));
CI95 = tinv([0.025 0.975], length(dQ3)-1);                    
dQ3CI95 = bsxfun(@times, sem, CI95(:)); 


sem = std(dQ4)/sqrt(length(dQ4));
CI95 = tinv([0.025 0.975], length(dQ4)-1);                    
dQ4CI95 = bsxfun(@times, sem, CI95(:)); 


figure;
errorbar(1,mean(dQ),dQCI95(1),dQCI95(2),'k')
hold on
errorbar(2,mean(dQ2),dQ2CI95(1),dQ2CI95(2),'k')
errorbar(3,mean(dQ3),dQ3CI95(1),dQ3CI95(2),'k')
errorbar(4,mean(dQ4),dQ4CI95(1),dQ4CI95(2),'k')
bar(1, mean(dQ),'facecolor',[0.5 0.7 0.8])
bar(2, mean(dQ2),'facecolor', [0.5 0.7 0.8])
bar(3, mean(dQ3),'facecolor', [0.5 0.7 0.8])
bar(4, mean(dQ4),'facecolor', [0.5 0.7 0.8])
xticks(1:4)
xticklabels({'Theta/Alpha' 'Beta' 'LowGamma' 'HighGamma'})
title(['Euclidean Distance from Center'])
set(gcf, 'color', 'w')
box off

[p,~,stats] = anova1(dQall, fIDq);
[p,~,stats] = anova1(diagQall, fIDq);

[p,~,stats] = anova1(dPall, fIDp);
[p,~,stats] = anova1(diagPall, fIDp);



%%
if ploton
    
   
    
    for x=1:27
        for j=1:4
            numInfluentialNodesQ(x,j) = sum(ssIDq==x&fIDq==j);
            numInfluentialNodesP(x,j) = sum(ssIDp==x&fIDp==j);
            
        end
    end
    figure;
    boxplot(numInfluentialNodesQ)
    title(['Number of Feature Nodes Qexp, nodeType ' num2str(nodeType)])
    xticklabels({'Theta/Alpha' 'Beta' 'LowGamma' 'HighGamma'})
    figure
    boxplot(numInfluentialNodesP)
    title(['Number of Feature Nodes POW, nodeType ' num2str(nodeType)])
    xticklabels({'Theta/Alpha' 'Beta' 'LowGamma' 'HighGamma'})
    
    meanQ = mean(numInfluentialNodesQ)
    medianQ = median(numInfluentialNodesQ)
    sumQ = sum(numInfluentialNodesQ)
    meanP = mean(numInfluentialNodesP)
    medianP = median(numInfluentialNodesP)
    sumP = sum(numInfluentialNodesP)
    table(numInfluentialNodesQ)
    table(numInfluentialNodesP)
    
    figure
    wc = wordcloud(gchQ);
    disp(['Nodal Q wordcloud top 5 count: ' num2str(wc.SizeData(1:5))])
    title(['Wordcloud Qexp Channels, nodeType ' num2str(nodeType)])
    
    figure
    wc1 = wordcloud(gchP);
    disp(['Nodal P wordcloud top 5 count: ' num2str(wc1.SizeData(1:5))])
    title(['Wordcloud POW Channels, nodeType ' num2str(nodeType)])
    
    figure;
    wc2 = wordcloud(mniQroi);
    disp(['MNI Q wordcloud top 5 count: ' num2str(wc2.SizeData(1:5))])
    title(['Wordcloud MNI ROI Qexp, nodeType ' num2str(nodeType)])
    
    figure;
    wc3 = wordcloud(mniProi);
    disp(['MNI P wordcloud top 5 count: ' num2str(wc3.SizeData(1:5))])
    title(['Wordcloud MNI ROI POW, nodeType ' num2str(nodeType)])
    
    figure;
    wc4 = wordcloud(yeoQroi);
    disp(['YEO Q wordcloud top 5 count: ' num2str(wc4.SizeData(1:5))])
    title(['Wordcloud YEO ROI Qexp, nodeType ' num2str(nodeType)])
    
    figure;
    wc5 = wordcloud(yeoProi);
    disp(['YEO P wordcloud top 5 count: ' num2str(wc5.SizeData(1:5))])
    title(['Wordcloud YEO ROI POW, nodeType ' num2str(nodeType)])
    
%     figure;
%     wc6 = wordcloud(agmniQ);
%     disp(['All Good Ch MNI Q wordcloud top 5 count: ' num2str(wc6.SizeData(1:5))])
%     
%     figure;
%     wc7 = wordcloud(agmniP);
%     disp(['All Good Ch MNI P wordcloud top 5 count: ' num2str(wc7.SizeData(1:5))])
%     
%     figure;
%     wc8 = wordcloud(agyeoQ);
%     disp(['All Good Ch YEO Q wordcloud top 5 count: ' num2str(wc8.SizeData(1:5))])
%     
%     figure;
%     wc9 = wordcloud(agyeoP);
%     disp(['All Good Ch YEO P wordcloud top 5 count: ' num2str(wc9.SizeData(1:5))])
    
    
end


% relWMlabelsQ = sum(strcmp('White_Matter',mniQroi))/sum(strcmp('White_Matter',agmniQ))*100;
% relNAlabelsQ = sum(strcmp('n/a',mniQroi))/sum(strcmp('n/a',agmniQ))*100;
% relGMlabelsQ = (sum(~strcmp('White_Matter',mniQroi)) - sum(strcmp('n/a',mniQroi)))/(sum(~strcmp('White_Matter',agmniQ)) - sum(strcmp('n/a',agmniQ)))*100;
% 
% 
% relWMlabelsP = sum(strcmp('White_Matter',mniProi))/sum(strcmp('White_Matter',agmniP))*100;
% relNAlabelsP = sum(strcmp('n/a',mniProi))/sum(strcmp('n/a',agmniP))*100;
% relGMlabelsP = (sum(~strcmp('White_Matter',mniProi)) - sum(strcmp('n/a',mniProi)))/(sum(~strcmp('White_Matter',agmniP)) - sum(strcmp('n/a',agmniP)))*100;

% figure;
% bar([relWMlabelsQ relWMlabelsP relGMlabelsQ relGMlabelsP relNAlabelsQ relNAlabelsP])
% xticklabels({'WM Q' 'WM P' 'GM Q' 'GM P' 'NA Q' 'NA P'})
% title(['Relative WM/GM/NA labels, nodeType ' num2str(nodeType)])

if ploton
figure;
[nQ, histQ] = cellhist(mniQroi);
title(['Influential Qexp Node ROI Counts, nodeType ' num2str(nodeType)])
figure;
[nP, histP] = cellhist(mniProi);
title(['Influential POW Node ROI Counts, nodeType ' num2str(nodeType)])
figure;
[nQa, histQa] = cellhist(agmniQ(ismember(agmniQ,mniQroi)));
title(['Node ROI Total Counts, nodeType ' num2str(nodeType)])
% figure;
% [nPa, histPa] = cellhist(agmniP(ismember(agmniP,mniProi)));
% title(['POW Node ROI Total Counts, nodeType ' num2str(nodeType)])
figure;
[nQf, histQf] = cellhist(mniQroi);
title(['Influential Qexp Node ROI Counts, nodeType ' num2str(nodeType)])

[ynQ, yhistQ] = cellhist(yeoQroi);
title(['Influential Qexp Node YEO Counts, nodeType ' num2str(nodeType)])
figure;
[ynP, yhistP] = cellhist(yeoProi);
title(['Influential POW Node YEO Counts, nodeType ' num2str(nodeType)])
figure;
[ynQa, yhistQa] = cellhist(agyeoQ(ismember(agyeoQ,yeoQroi)));
title(['Node YEO Total Counts, nodeType ' num2str(nodeType)])
figure;
% [ynPa, yhistPa] = cellhist(agyeoP(ismember(agyeoP,yeoProi)));
% title(['POW Node YEO Total Counts, nodeType ' num2str(nodeType)])



for i=1:10
relPercQ(i) = sum(strcmp(histQ(i),mniQroi))/sum(strcmp(histQ(i),agmniQ))*100;
end
figure;
[~,si] = sort(relPercQ,'descend');
bar(relPercQ(si))
chistQ = histQ(1:10);
xticks(1:10)
xticklabels(chistQ(si))
xtickangle(45)
title(['Relative Frequency Influential Contribution Qexp, nodeType ' num2str(nodeType)])

for i=1:10
relPercP(i) = sum(strcmp(histP(i),mniProi))/sum(strcmp(histP(i),agmniP))*100;
end
figure;
[~,si] = sort(relPercP,'descend');
chistP = histP(1:10);
bar(relPercP(si))
xticks(1:10)
xticklabels(chistP(si))
xtickangle(45)
title(['Relative Frequency Influential Contribution POW, nodeType ' num2str(nodeType)])

for i=1:length(yhistQ)
yrelPercQ(i) = sum(strcmp(yhistQ(i),yeoQroi))/sum(strcmp(yhistQ(i),agyeoQ))*100;
end
figure;
[~,si] = sort(yrelPercQ,'descend');
bar(yrelPercQ(si))
chistQ = yhistQ(1:length(yhistQ));
xticks(1:length(yhistQ))
xticklabels(chistQ(si))
xtickangle(45)
title(['Relative Frequency Influential Yeo Network Contribution Qexp, nodeType ' num2str(nodeType)])

for i=1:length(yhistP)
yrelPercP(i) = sum(strcmp(yhistP(i),yeoProi))/sum(strcmp(yhistP(i),agyeoP))*100;
end
figure;
[~,si] = sort(yrelPercP,'descend');
chistP = yhistP(1:length(yhistP));
bar(yrelPercP(si))
xticks(1:length(yhistP))
xticklabels(chistP(si))
xtickangle(45)
title(['Relative Frequency Influential Yeo Network Contribution POW, nodeType ' num2str(nodeType)])

end