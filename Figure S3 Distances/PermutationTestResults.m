clear variables
close all

% path='\\lnlogin04.res.utc.com\hpwip_tmp2_usp\systems\project\pw_blade_inspection\AdditiveManufacturing\';
% inputFolder='dissimilarityMeasures/dataUse/';
% outputFolder='dissimilarityMeasures/testFullUpdated/';
chrNo=0; % this is for full genome
nHist=10;

% mesIdUse=[1 3 5]; % if you just want to plot spectral distance, change  mesIdUse=[3]
mesIdUse=[3]
grIdUse=[ 5];

disNames{1}='D_H';
disNames{end+1}='D_J';
disNames{end+1}='D_{\lambda}';
disNames{end+1}='D_{\tau}';
disNames{end+1}='D_{\Delta}';
disNamesUse=disNames(mesIdUse);

%% read
if 1
    load randomPertGraphs dnullFull  res1 % this is full genome
    graphRepresentTypeUse=res1.graphRepresentTypeUse(grIdUse);
    measureTypeUse=res1.measureTypeUse(mesIdUse);
    titleTypeUse=res1.titleTypeUse(grIdUse);
else
    outP=[path outputFolder];
    fnames=dir([outP 'measure*']);
    res1=load([outP fnames(1).name],'measureTypeUse','graphRepresentTypeUse','titleType');
    graphRepresentTypeUse=res1.graphRepresentTypeUse(grIdUse);
    measureTypeUse=res1.measureTypeUse(mesIdUse);
    titleTypeUse=res1.titleType(grIdUse);
    
    for j=1:length(fnames)
        j
        res=load([outP fnames(j).name],'dnull');
        dnullFull(:,:,j)=res.dnull;
    end
end

%%
gmfib=load('fibGMPointComparision_CombinedFib1234','dfibgm','measureTypeUse','graphRepresentTypeUse');
graphR=gmfib.graphRepresentTypeUse;
measureT=gmfib.measureTypeUse;
%dfibgm: crNo x measureType x graph type

indGR=[];
for i=1:length(graphRepresentTypeUse)
    for j=1:length(graphR)
        if strcmp(graphRepresentTypeUse{i},graphR{j})
            indGR(i)=j;
        end
    end
end

indMes=[];
for i=1:length(measureTypeUse)
    for j=1:length(measureT)
        if strcmp(measureTypeUse{i},measureT{j})
            indMes(i)=j;
        end
    end
end
dfibgm=gmfib.dfibgm(chrNo+1,indMes,indGR);
dfibgm=reshape(dfibgm,size(dfibgm,2),size(dfibgm,3));

%%
% for k=1:length(graphRepresentTypeUse)
%     gSbaseFib1=createGraphRepresentation(baseHGFib,graphRepresentTypeUse{k});
%     gSTestGM=createGraphRepresentation(testHGGM,graphRepresentTypeUse{k});  
%     for i=1:length(measureTypeUse)
%         dfibgm(i,k)=graphDisMeasures(gSbaseFib1,gSTestGM,measureTypeUse{i},params);
%     end  
% end

for i=1:length(measureTypeUse)
    if strcmp(measureTypeUse{i},'lamdaL2')
        measureTypeUse{i}='Spectral';
    end
end
%%
figure
for k=1:length(graphRepresentTypeUse)
    for i=1:length(measureTypeUse)
        dnullUse=dnullFull(mesIdUse(i),grIdUse(k),:);       
        pfibgm(i,k)=sum(dnullUse(:)>dfibgm(i,k))/length(dnullUse(:));
        subplot(length(graphRepresentTypeUse),length(measureTypeUse),i)
        hist(dnullUse(:),nHist)
        hold on
        %lLeg(1)=plot(dfibgm(i,k),0,'*');
        ttl=title(measureTypeUse{i})
        xlb=xlabel(disNamesUse{i})
        if i==1
            ylb=ylabel('Frequency')
            %legend(lLeg,{'Fib-GM'},'Location','Best')
        end
        str=sprintf('p=%1.3f',pfibgm(i,k));
        %title(str)
        axis tight
        set(ttl,'fontsize',16)
        set(ylb,'fontsize',16)
        set(xlb,'fontsize',16)
        ax=gca;
        ax.LineWidth=2;
        ax.Box='on';
        ax.FontSize=12;
        %ax.XTickLabel=[];
        %ax.YTickLabel=[];
    end
end

%% Stem plot - SL
figure
[pdf_edge_empirical, boundary_edge]= histcounts(dnullUse(:),30); 
pts_edge = ( boundary_edge(1:(end-1)) + boundary_edge(2:end) )/2;
stem(pts_edge,pdf_edge_empirical,'-b','LineWidth',1.5,'Marker','o','MarkerSize',5); hold on;
stem(dfibgm,0.5*max(pdf_edge_empirical), 'LineStyle','--','Color','r','LineWidth',2); 
if dfibgm < .01
    text(dfibgm,0.55*max(pdf_edge_empirical),'Pvalue < 0.01','HorizontalAlignment', 'right');
else
    text(dfibgm,0.55*max(pdf_edge_empirical),...
        sprintf('Pvalue = %1.2f',floor(dfibgm*100)/100),'HorizontalAlignment', 'right');
end
set(gca,'YLim',get(gca,'YLim')*1.1,'TickLength',[0 0])
% xlim([7.95E-4 8.25E-4])
% title('Randomly sampled genes')
