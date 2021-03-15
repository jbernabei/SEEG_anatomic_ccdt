% CCDT network state space
% Per subject, network state space is the real-time combination of Qexp and
% Pow values across feature selected nodes across frequency bands. 

% Run after CCDTnodeAnalyze2.m


% VB 01/21
netSSnodes=cell(27,1);
for i = 1:27
    for f = 1:4
        
        tQ{f} = gchQ(ssIDq==i&fIDq==f);
        tP{f} = gchP(ssIDp==i&fIDp==f);
        
    end
    
    netSSnodes{i} = [tQ tP];
end


netSS = cell(27,2);
for i=1:27
    for f=1:4
        chIndQ = netSSnodes{i}{f};
        chIndP = netSSnodes{i}{f+4};
        ssChQ = ismember(NFstruct(i).gchlbl,chIndQ);
        ssChP = ismember(NFstructP(i).gchlbl(:,1),chIndP);
        
        czQ = zscore(NFstruct(i).fbands(f).NFqexp_nodal,1);
        czP = zscore(NFstructP(i).fbands(f).NFpow,1);
        
        cLT = LTqexp{i}(ssChQ,f,1);
        cLTP = LTpow{i}(ssChP,f,1);

        tQss{f} = czQ(ssChQ,:);
        tPss{f} = czP(ssChP,:);
        
        tQb{f} = cLT;
        tPb{f} = cLTP;
    end

    tnetSS = [tQss tPss];
    tBval = [tQb tPb];
    for j=1:8
    netSS{i,1} = [netSS{i,1} tnetSS{j}'];
    netSS{i,2} = [netSS{i,2} tBval{j}'];
    end
end



%% plots

for i=[14,16]
    figure(i);clf
    cvRT = behavStruct(i).vRT
    [a,b] = sort(cvRT, 'descend');
    cBval = netSS{i,2};
    [~,c] = sort(cBval, 'ascend');
    subplot(3,1,1)
    imagesc(netSS{i}(b(a>0&a<999),c)')
    ylabel('Feature Node (sorted B val)')
    caxis([-2 2])
%     subplot(4,1,2)
%     [~, maxi] = max(netSS{i}(b(a>0&a<999),c)');
%     plot(maxi.*-1,'.')
%     lsline
%     box off
%     axis tight
    subplot(3,1,3)
    plot(cvRT(b(a>0&a<999)),'k.','markersize',12)
    ylabel('RT')
    xlabel('Trial (sorted RT)')
    axis tight
    box off
    set(gcf,'color','w')
    print(['netSS.ps'],'-dpsc2','-bestfit','-append');
    
    feature_space = netSS{i}(b(a>0&a<999),c);
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(feature_space);
    subplot(3,1,2)
    [RHO,PVAL] = corr([1:length(SCORE(:,2))]',SCORE(:,2));
    
    plot([1:length(SCORE(:,2))]',SCORE(:,2),'k.')
    title(sprintf('Pearson corr p = %2f',PVAL))
    ylabel('First principal component score')
    xlabel('Trial (sorted RT)')
    lsline
    
end


