% Figures and Analyses for BlindnessNOS structure paper

clear all;
close all;

% This is a hard-coded list of indexes for the different blind and sighted
% populations in the data. It is critical that these align with the
% measures.

indexsight=linspace(60,118,59);
indexblind=linspace(1,59,59);
indexanophthalmic=linspace(1,6,6);
indexcongenital=linspace(7,27,21);
indexpostnatal=linspace(28,40,13);
indexlca=linspace(41,59,19);

indexpostnatal_beforepuberty= indexpostnatal([1 3 6 7 10 12 13]);
indexpostnatal_afterpuberty =indexpostnatal([2 4 5 8 9 11]);

indexlca1=linspace(41,46,6);
indexrpe65=linspace(47,52,6);
indexcrb1=linspace(53,57,5);
indexcep290=linspace(58,59,2);



% The test / retest data are for the visual data set only

indextest=linspace(1,24,24);
indexretest=linspace(25,48,24);

% This is the hard-coded threshold for excluding outliers, in units of
%  standard deviations.

OutThresh=4.5;


% Create a header to mark out the text output

fprintf('\n\n\n***************************************************\n\n')

% load the data. This includes a matrix of measures for the visual and
%  auditory pathway.

load('/Users/Shared/Matlab/Aguirre Lab Repo/BlindnessNOS_data.mat');

%%%%%%%%%%%%%%%%%%%
% Change these lines to flip the analysis from visual to auditory data
%

fprintf('\n\nVISUAL DATA ANALYSIS\n\n');
data=visual_data;
%OutFileStem='/Users/aguirre/Temp/Structure_Function_Alterations_Blinds/Figures_and_Tables_v4.0/BlindnessNOSFigs_Visual/';
OutFileStem='/Users/aguirre/Temp/';


%fprintf('\n\nAUDITORY DATA ANALYSIS\n\n');
%data=auditory_data;
%OutFileStem='/Users/aguirre/Dropbox/Structure_Function_Alterations_Blinds/Figures_and_Tables_v4.0/BlindnessNOSFigs_Auditory/';

%
%
%%%%%%%%%%%%%%%%%%%


sizer=size(data);
NumSubjects=sizer(1);
NumMeasures=sizer(2);

% Place to save figures

FigIndex=1;

% Create a little subject demographics table

fprintf('\n\nSubject Demographics:\n\n');
fprintf(['Sighted: ' num2str(sum(gender(indexsight))) 'm / ' num2str(length(indexsight)-sum(gender(indexsight))) 'f, age= ' num2str(mean(ages(indexsight)),'%1.0f') ' ± ' num2str(std(ages(indexsight)),'%1.0f') '\n']);
fprintf(['All blind: ' num2str(sum(gender(indexblind))) 'm / ' num2str(length(indexblind)-sum(gender(indexblind))) 'f, age= ' num2str(mean(ages(indexblind)),'%1.0f') ' ± ' num2str(std(ages(indexblind)),'%1.0f') '\n']);
fprintf(['  Anophthalmic: ' num2str(sum(gender(indexanophthalmic))) 'm / ' num2str(length(indexanophthalmic)-sum(gender(indexanophthalmic))) 'f, age= ' num2str(mean(ages(indexanophthalmic)),'%1.0f') ' ± ' num2str(std(ages(indexanophthalmic)),'%1.0f') '\n']);
fprintf(['  Other early: ' num2str(sum(gender(indexcongenital))) 'm / ' num2str(length(indexcongenital)-sum(gender(indexcongenital))) 'f, age= ' num2str(mean(ages(indexcongenital)),'%1.0f') ' ± ' num2str(std(ages(indexcongenital)),'%1.0f') '\n']);
fprintf(['  Leber amaurosis: ' num2str(sum(gender(indexlca))) 'm / ' num2str(length(indexlca)-sum(gender(indexlca))) 'f, age= ' num2str(mean(ages(indexlca)),'%1.0f') ' ± ' num2str(std(ages(indexlca)),'%1.0f') '\n']);
fprintf(['  Other late: ' num2str(sum(gender(indexpostnatal))) 'm / ' num2str(length(indexpostnatal)-sum(gender(indexpostnatal))) 'f, age= ' num2str(mean(ages(indexpostnatal)),'%1.0f') ' ± ' num2str(std(ages(indexpostnatal)),'%1.0f') '\n']);
fprintf('\n\n');

% Report age of onset of blindness for postnatal and congential:

fprintf('\n\nMean age of blindness onset ± SD [range]:\n\n');
fprintf(['  Other early: ' num2str(mean(blindness_onset_max_years(indexcongenital,1))) '  [' num2str(min(blindness_onset_max_years(indexcongenital,1))) ' - ' num2str(max(blindness_onset_max_years(indexcongenital,1))) ']\n']);
fprintf(['  Other late: ' num2str(mean(blindness_onset_max_years(indexpostnatal,1)),'%1.0f') ' ± '  num2str(std(blindness_onset_max_years(indexpostnatal,1)),'%1.0f')  '  [' num2str(min(blindness_onset_max_years(indexpostnatal,1))) ' - ' num2str(max(blindness_onset_max_years(indexpostnatal,1))) ']\n']);
fprintf('\n\n');


% Report age at max blindness for postnatal and congential:

fprintf('\n\nMean age of max blindness ± SD [range]:\n\n');
fprintf(['  Other early: ' num2str(mean(blindness_onset_max_years(indexcongenital,2)),'%1.0f') '  [' num2str(min(blindness_onset_max_years(indexcongenital,2))) ' - ' num2str(max(blindness_onset_max_years(indexcongenital,2))) ']\n']);
fprintf(['  Other late: ' num2str(mean(blindness_onset_max_years(indexpostnatal,2)),'%1.0f') ' ± '  num2str(std(blindness_onset_max_years(indexpostnatal,2)),'%1.0f')  '  [' num2str(min(blindness_onset_max_years(indexpostnatal,2))) ' - ' num2str(max(blindness_onset_max_years(indexpostnatal,1))) ']\n']);
fprintf('\n\n');


% Number of years of max blindness:

fprintf('\n\nNumber of years of max blindness ± SD [range]:\n\n');
fprintf(['  Other early: ' num2str(mean(ages(indexcongenital)-blindness_onset_max_years(indexcongenital,2)),'%1.0f') '  [' num2str(min(ages(indexcongenital)-blindness_onset_max_years(indexcongenital,2))) ' - ' num2str(max(ages(indexcongenital)-blindness_onset_max_years(indexcongenital,2))) ']\n']);
fprintf(['  Other late: ' num2str(mean(ages(indexpostnatal)-blindness_onset_max_years(indexpostnatal,2)),'%1.0f') ' ± '  num2str(std(ages(indexpostnatal)-blindness_onset_max_years(indexpostnatal,2)),'%1.0f')  '  [' num2str(min(ages(indexpostnatal)-blindness_onset_max_years(indexpostnatal,2))) ' - ' num2str(max(ages(indexpostnatal)-blindness_onset_max_years(indexpostnatal,1))) ']\n']);
fprintf('\n\n');



% Report the mean and SEM of the supratentorial volume by group

fprintf('\n\nMean (±SEM) supertentorial volume [liters]:\n\n');
fprintf(['Sighted: ' num2str(mean(super_tentorial_volume(indexsight))*1e-6,'%1.2f') ' ± ' num2str(std(super_tentorial_volume(indexsight))/sqrt(length(indexsight))*1e-6,'%1.2f') '\n']);
fprintf(['All blind: ' num2str(mean(super_tentorial_volume(indexblind))*1e-6,'%1.2f') ' ± ' num2str(std(super_tentorial_volume(indexblind))/sqrt(length(indexblind))*1e-6,'%1.2f') '\n']);
fprintf(['  Anophthalmic: ' num2str(mean(super_tentorial_volume(indexanophthalmic))*1e-6,'%1.2f') ' ± ' num2str(std(super_tentorial_volume(indexanophthalmic))/sqrt(length(indexanophthalmic))*1e-6,'%1.2f') '\n']);
fprintf(['  Other early: ' num2str(mean(super_tentorial_volume(indexcongenital))*1e-6,'%1.2f') ' ± ' num2str(std(super_tentorial_volume(indexcongenital))/sqrt(length(indexcongenital))*1e-6,'%1.2f') '\n']);
fprintf(['  Leber amaurosis: ' num2str(mean(super_tentorial_volume(indexlca))*1e-6,'%1.2f') ' ± ' num2str(std(super_tentorial_volume(indexlca))/sqrt(length(indexlca))*1e-6,'%1.2f') '\n']);
fprintf(['  Other late: ' num2str(mean(super_tentorial_volume(indexpostnatal))*1e-6,'%1.2f') ' ± ' num2str(std(super_tentorial_volume(indexpostnatal))/sqrt(length(indexpostnatal))*1e-6,'%1.2f') '\n']);
fprintf('\n\n');

% Report the number of subjects who have each kind of "other" measure

fprintf('\n\nNumber of subjects with additional brain measures:\n\n');
fprintf( '         CBF  Sentences FA\n');
fprintf(['Sighted: ' num2str(length(indexsight)-sum(isnan(othermeasures(indexsight,:)))) '\n']);
fprintf(['Blind:   ' num2str(length(indexblind)-sum(isnan(othermeasures(indexblind,:)))) '\n']);
fprintf('\n\n');




% Report the Cohen's d effect size of blind vs sighted for each measure

fprintf('Cohen-s d effect sizes for blind vs. sighted for each measure:\n');
for i=1:NumMeasures
    d=(mean(data(indexsight,i))-mean(data(indexblind,i)))/ ( ( std(data(indexsight,i)) + std(data(indexblind,i)) ) /2);
    fprintf(['measure ' num2str(i) ': ' num2str(d,'%1.1f') '\n']);
end
fprintf('\n\n');

%  fprintf('Number of SEMs between the blind vs. sighted for each measure:\n');
%  for i=1:NumMeasures
%      d=(mean(data(indexsight,i))-mean(data(indexblind,i)))/ ( ( std(data(indexsight,i))/sqrt(length(indexsight)) + std(data(indexblind,i))/sqrt(length(indexblind)) ) /2);
%      fprintf(['measure ' num2str(i) ': ' num2str(d) '\n']);
%  end
%  fprintf('\n\n');





% Create parametric box plots -- limited range

figtmp=figure('name','Sighted Scores');
notBoxPlot(data(indexsight,:));
ylim([-3 3]);
title('Sighted scores, mean, ±1.96 SEM, ±1SD');
pbaspect([2 2 1])
box off;
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;

figtmp=figure('name','Blind Scores');
notBoxPlot(data(indexblind,:));
ylim([-3 3]);
title('Blind scores, mean, ±1.96 SEM, ±1SD');
pbaspect([2 2 1])
box off;
hold off
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


% Create parametric box plots -- full range

figtmp=figure('name','Sighted Scores');
notBoxPlot(data(indexsight,:));
ylim([-5 5]);
title('Sighted scores, mean, ±1.96 SEM, ±1SD');
pbaspect([2 2 1])
box off;
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;

figtmp=figure('name','Blind Scores');
notBoxPlot(data(indexblind,:));
ylim([-5 5]);
title('Blind scores, mean, ±1.96 SEM, ±1SD');
pbaspect([2 2 1])
box off;
hold off
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;





% Test if the standard deviations differ between the two groups for each
% measure
% Commented out as this is not used in the paper

% fprintf('Test if the standard deviations differ between the two groups for each measure:\n');
%
% for i=1:NumMeasures
%     [p,stats]=vartestn([data(indexsight,i),data(indexblind,i)],'Display','off');
%     fprintf(['measure dim ',num2str(i),'- c2(1, N = ' num2str(length(indexsight)+length(indexblind)) ') = ',num2str(stats.chisqstat),', p=',num2str(p)]);
%     fprintf('\n');
% end
% fprintf('\n\n');




% Create the cross-correlation plots

M=corr(data(indexsight,:));
for i=1:NumMeasures M(i,i)=0; % set the diagonal to zero, as opposed to one
end
figtmp=figure('name','Sighted measure correlations');
h = imagesc(M);
axis square;
caxis([-0.6 0.6]);
C=polarmap();
colormap(C);
colorbar;
set(gca, 'visible', 'off') ;
title('Sighted measure correlations');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;

M=corr(data(indexblind,:));
for i=1:NumMeasures M(i,i)=0; % set the diagonal to zero, as opposed to one
end
figtmp=figure('name','Blind measure correlations');
h = imagesc(M);
axis square;
caxis([-0.6 0.6]);
C=polarmap();
colormap(C);
colorbar;
set(gca, 'visible', 'off') ;
title('Blind measure correlations');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;

M=corr(data);
for i=1:NumMeasures M(i,i)=0; % set the diagonal to zero, as opposed to one
end
figtmp=figure('name','All subjects measure correlations');
h = imagesc(M);
axis square;
caxis([-0.6 0.6]);
C=polarmap();
colormap(C);
colorbar;
set(gca, 'visible', 'off') ;
title('All subjects (blinda and sighted) measure correlations');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;




%%%%%%%%%%%%%%%%%%%%
% Run some analyses to determine the amount of variance each data set
%  can explain in other datasets
%%%%%%%%%%%%%%%%%%%%%%


ExplainedDataLabels={'define sight -> test sight' 'define sight -> test blind' 'define sight -> test all';...
    'define blind -> test sight' 'define blind -> test blind' 'define blind -> test all';...
    'define all -> test sight' 'define all -> test blind' 'define all -> test all'};
ExplainedData=zeros(3,3);





% conduct the PCA analysis with the sighted subjects only

[sight_coeff,temp_score,~,~,~,~] = loo_pca(data(indexsight,:));

% Calculate the PCA scores for the blind (and again for the sighted) using
% the sighted coeff

for i=1:NumSubjects
    for component=1:NumMeasures
        sight_score(i,component)=sum(sight_coeff(:,component)'.*data(i,:));
    end
end

ExplainSight=0;
ExplainBlind=0;
ExplainAll=0;

for i=1:3
    ExplainSight=ExplainSight+(std(sight_score(indexsight,i))^2)/(sum(std(sight_score(indexsight,:)).^2));
    ExplainBlind=ExplainBlind+(std(sight_score(indexblind,i))^2)/(sum(std(sight_score(indexblind,:)).^2));
    ExplainAll=ExplainAll+(std(sight_score(:,i))^2)/(sum(std(sight_score(:,:)).^2));
end

ExplainedData(1,1)=ExplainSight;
ExplainedData(1,2)=ExplainBlind;
ExplainedData(1,3)=ExplainAll;





% conduct the PCA analysis with the blind subjects only

[blind_coeff,~,~,~,~,~] = loo_pca(data(indexblind,:));

% Calculate the PCA scores for the sighted (and again for the blind) using
% the blind coeffs

for i=1:NumSubjects
    for component=1:NumMeasures
        blind_score(i,component)=sum(blind_coeff(:,component)'.*data(i,:));
    end
end

ExplainSight=0;
ExplainBlind=0;
ExplainAll=0;

for i=1:3
    ExplainSight=ExplainSight+(std(blind_score(indexsight,i))^2)/(sum(std(blind_score(indexsight,:)).^2));
    ExplainBlind=ExplainBlind+(std(blind_score(indexblind,i))^2)/(sum(std(blind_score(indexblind,:)).^2));
    ExplainAll=ExplainAll+(std(blind_score(:,i))^2)/(sum(std(blind_score(:,:)).^2));
end

ExplainedData(2,1)=ExplainSight;
ExplainedData(2,2)=ExplainBlind;
ExplainedData(2,3)=ExplainAll;






% conduct the PCA analysis with all subjects

[all_coeff,~,~,~,~,~] = loo_pca(data);

% Calculate the PCA scores for the sighted (and again for the blind) using
% the all coeffs

for i=1:NumSubjects
    for component=1:NumMeasures
        all_score(i,component)=sum(all_coeff(:,component)'.*data(i,:));
    end
end

ExplainSight=0;
ExplainBlind=0;
ExplainAll=0;

for i=1:3
    ExplainSight=ExplainSight+(std(all_score(indexsight,i))^2)/(sum(std(all_score(indexsight,:)).^2));
    ExplainBlind=ExplainBlind+(std(all_score(indexblind,i))^2)/(sum(std(all_score(indexblind,:)).^2));
    ExplainAll=ExplainAll+(std(all_score(:,i))^2)/(sum(std(all_score(:,:)).^2));
end

ExplainedData(3,1)=ExplainSight;
ExplainedData(3,2)=ExplainBlind;
ExplainedData(3,3)=ExplainAll;


% Make a bar plot reporting these results


figtmp=figure('name','Proportion of variance explained by 3 PCs defined in different groups');

bar(ExplainedData)
ylim([0,1])
set(gca,'XTickLabel',ExplainedDataLabels)
hold on
plot(xlim,[0.333 0.3333], 'r')
pbaspect([2 2 1])
title('Proportion of variance explained by 3 PCs defined in different groups');
box off;
hold off

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;

% Dump out the numbers

fprintf('The proportion of variance explained in each group by PC1-3 defined\n');
fprintf('  using that group or defined using another group. These are the\n');
fprintf('  values in the bar plot.\n');
ExplainedData
fprintf('\n\n');






% Specifically report how much variance the mean difference between blind
%  and sighted explains in the sighted and blind groups alone

fprintf('The proportion of variance explained in each group by the mean\n');
fprintf('  difference between blind and sighted groups:\n');

GroupMeanDiffs=mean(data(indexsight,:),1)-mean(data(indexblind,:),1);

for i=1:NumSubjects
    GroupMeanDiff_score(i)=sum(GroupMeanDiffs.*data(i,:));
end

GroupMeanExplainSight=0;
GroupMeanExplainBlind=0;
GroupMeanExplainAll=0;

for i=1:NumMeasures
    GroupMeanExplainSight=GroupMeanExplainSight+nancorr(GroupMeanDiff_score(indexsight),data(indexsight,i))^2;
    GroupMeanExplainBlind=GroupMeanExplainBlind+nancorr(GroupMeanDiff_score(indexblind),data(indexblind,i))^2;
    GroupMeanExplainAll=GroupMeanExplainAll+nancorr(GroupMeanDiff_score(:),data(:,i))^2;
end

GroupMeanExplainSight=GroupMeanExplainSight/sum(std(data,1).^2);
GroupMeanExplainBlind=GroupMeanExplainBlind/sum(std(data,1).^2);
GroupMeanExplainAll=GroupMeanExplainAll/sum(std(data,1).^2);

fprintf(['Explain sighted: ' num2str(GroupMeanExplainSight,'%1.2f') '\n']);
fprintf(['Explain blind: ' num2str(GroupMeanExplainBlind,'%1.2f') '\n']);
fprintf(['Explain all: ' num2str(GroupMeanExplainAll,'%1.2f') '\n']);
fprintf('\n\n');

% Plot these data

figtmp=figure('name','Explained variance by PC1 of Both');
titles={'sighted' 'blind' 'all'};

subplot(3,1,1);
bar(3,GroupMeanExplainSight*100,'k');
hold on
xlim([1,NumMeasures]);
ylim([0,50]);
pbaspect([2 2 1])
title([titles{1}]);
box off;
hold off

subplot(3,1,2);
bar(3,GroupMeanExplainBlind*100,'r');
hold on
xlim([1,NumMeasures]);
ylim([0,50]);
pbaspect([2 2 1])
title([titles{1}]);
box off;
hold off

subplot(3,1,3);
bar(3,GroupMeanExplainAll*100,'k');
hold on
xlim([1,NumMeasures]);
ylim([0,50]);
pbaspect([2 2 1])
title([titles{1}]);
box off;
hold off

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;









% Report how well the sighted variation accounts for the effects of
% blindness

fprintf('Regress the mean sight-blind differences on the first two PCs\n');
fprintf('  derived from the sighted:\n');

[~,~,~,~,stats] = regress(GroupMeanDiffs',[sight_coeff(:,1:2) ones(9,1)]);

stats

1-fcdf(stats(2),2,6)

fprintf('\n\n');

fprintf('Regress the mean sight-blind differences on the first three PCs\n');
fprintf('  derived from the sighted:\n');

[~,~,~,~,stats] = regress(GroupMeanDiffs',[sight_coeff(:,1:3) ones(9,1)]);

stats

1-fcdf(stats(2),3,6)

fprintf('\n\n');





% Create a figure that gives the explained variance plot with SEM bars
%  from resampling for blind, sighted, and all data

[ ~, ~, mean_explained_sight, sem_explained_sight ] = boot_pca(data(indexsight,:));
[ ~, ~, mean_explained_blind, sem_explained_blind ] = boot_pca(data(indexblind,:));
[ ~, ~, mean_explained_all, sem_explained_all ] = boot_pca(data);

xvals=linspace(1,NumMeasures,NumMeasures);

figtmp=figure('name','Explained variance by PC ±SEM');
titles={'sighted' 'blind' 'all'};

subplot(3,1,1);
plot(xvals,mean_explained_sight,'.k');
hold on
[temp_fit]=fit(xvals',mean_explained_sight,'smoothingspline');
plot(temp_fit,'-k')
hold on
[temp_fit]=fit(xvals',mean_explained_sight+sem_explained_sight,'smoothingspline');
plot(temp_fit,'-b')
hold on
[temp_fit]=fit(xvals',mean_explained_sight-sem_explained_sight,'smoothingspline');
plot(temp_fit,'-b')
xlim([1,NumMeasures]);
ylim([0,50]);
pbaspect([2 2 1])
title([titles{1}]);
box off;
hold off

subplot(3,1,2);
plot(xvals,mean_explained_blind,'.r');
hold on
[temp_fit]=fit(xvals',mean_explained_blind,'smoothingspline');
plot(temp_fit,'-r')
hold on
[temp_fit]=fit(xvals',mean_explained_blind+sem_explained_blind,'smoothingspline');
plot(temp_fit,'-b')
hold on
[temp_fit]=fit(xvals',mean_explained_blind-sem_explained_blind,'smoothingspline');
plot(temp_fit,'-b')
xlim([1,NumMeasures]);
ylim([0,50]);
pbaspect([2 2 1])
title([titles{2}]);
box off;
hold off

subplot(3,1,3);
plot(xvals,mean_explained_all,'.k');
hold on
[temp_fit]=fit(xvals',mean_explained_all,'smoothingspline');
plot(temp_fit,'-k')
hold on
[temp_fit]=fit(xvals',mean_explained_all+sem_explained_all,'smoothingspline');
plot(temp_fit,'-b')
hold on
[temp_fit]=fit(xvals',mean_explained_all-sem_explained_all,'smoothingspline');
plot(temp_fit,'-b')
xlim([1,NumMeasures]);
ylim([0,50]);
pbaspect([2 2 1])
title([titles{3}]);
box off;
hold off

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;










% conduct the PCA analysis with all subjects. This is done as a LOO to
%  obtain the coefficients (and scores). The first three components are
%  then subjected to an orthogonal procrustes rotation to align the
%  coefficients with a target matrix.
%  The target matrix in this case is the first three axes defined by the
%  blind-only PCA analysis.

Target=pca(data(indexblind,:));

coeff=pca(data);

% Calculate the transformation matrix for the first three principal
% componets from the all subjects data to align with the blind axis

[~,T]= rotatefactors(coeff(:,1:3),'Method','procrustes','Type','orthogonal','Target',Target(:,1:3));
coeff(:,1:3)=coeff(:,1:3)*T;
coeff(:,4:end)=NaN;

% Flip and swap some dimensions to order properly for plotting

coeff(:,1)=coeff(:,1)*(-1);
coeff(:,[1,2])=coeff(:,[2,1]);


% Obtain the subject scores on the dimensions by LOO

[~,all_score,~,~,~,~] = loo_pca(data,T,3);


% Flip and swap some dimensions to order properly for plotting

all_score(:,1)=all_score(:,1)*(-1);
all_score(:,[1,2])=all_score(:,[2,1]);


% Now perform a boot strap analysis to obtain the SEM on the coefficients

[ boot_all_coeff, sem_all_coeff, ~, ~ ] = boot_pca(data, T,3);


% Swap some dimensions to order properly for plotting

sem_all_coeff(:,[1,2])=sem_all_coeff(:,[2,1]);





% Create a matrix plot of the loading on the significant anatomical
%  dimensions in the first three principal components

M=coeff;
noonsig=find(abs(coeff./sem_all_coeff)<1.95);
M(noonsig)=0;

figtmp=figure('name','Anatomical weights matrix, significant values only');
h = imagesc(M);
axis square;
caxis([-0.7 0.7]);
C=polarmap();
colormap(C);
colorbar;
set(gca, 'visible', 'off') ;
title('Anatomical weights on 3 dimensions');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;




% Create a matrix plot of the loading on all anatomical
%  dimensions in the first three principal components

M=coeff;

figtmp=figure('name','Anatomical weights matrix all values');
h = imagesc(M);
axis square;
caxis([-0.7 0.7]);
C=polarmap();
colormap(C);
colorbar;
set(gca, 'visible', 'off') ;
title('Anatomical weights on 3 dimensions');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;







% Calculate the test / re-test reliability for the 24 subjects
% with test / re-test anatomical measures. We will determine the
% test / re-test for each of the 3 studied dimensions.

for i=1:48
    for component=1:NumMeasures
        test_retest_score(i,component)=sum(coeff(:,component)'.*visual_test_retest(i,:));
    end
end

YRanges=[-2 4; -2 4; -2 4];

figtmp=figure('name','Test - Retest reliability for the three measurement dimensions');
titles={'Dim1' 'Dim2' 'Dim3'};

for i=1:3
    subplot(3,1,i);
    scatter(test_retest_score(indextest,i),test_retest_score(indexretest,i),40,[.5 0 .75],'+');
    hold on
    ylim(YRanges(i,:));
    xlim(YRanges(i,:));
    pbaspect([2 2 1])
    title([titles{i}]);
    box off;
    hold off
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;

% Report the test / re-test reliability of these measures

fprintf('In a population of 24 subjects, the test / re-test reliability of scores on \n');
fprintf('  each of the three anatomical dimensions was:\n');
for i=1:3
    fprintf(['Dim ' num2str(i) ': ' num2str(nancorr(test_retest_score(indextest,i),test_retest_score(indexretest,i)),'%1.2f') '\n']);
end
fprintf('\n\n');







% Make a little table giving the amount of variance each component in this
% space explains in the sighted, blind, and combined datasets

% Calculate the SEM by bootstrap resampling within group.

AllIndicesToSample=linspace(1,NumSubjects,NumSubjects);

nBootstraps = 1000;

BootstrapExplainedAll = zeros(3,nBootstraps);
BootstrapExplainedSight = zeros(3,nBootstraps);
BootstrapExplainedBlind = zeros(3,nBootstraps);

for b = 1:nBootstraps
    
    % First do a resample from all the subjects
    
    % This will sample with replacement from the rows of X,
    % matching the total number of available rows in X.
    
    SampleWithReplace= ...
        datasample(linspace(1,NumSubjects, ...
        NumSubjects),NumSubjects,'Replace',true);
    
    RandomSampleIndicesAll = AllIndicesToSample(SampleWithReplace);
    
    % Now the blind
    
    SampleWithReplace= ...
        datasample(linspace(1,length(indexblind), ...
        length(indexblind)),length(indexblind),'Replace',true);
    
    RandomSampleIndicesBlind = indexblind(SampleWithReplace);
    
    % Now the sighted
    
    SampleWithReplace= ...
        datasample(linspace(1,length(indexsight), ...
        length(indexsight)),length(indexsight),'Replace',true);
    
    RandomSampleIndicesSight = indexsight(SampleWithReplace);
    
    
    for i=1:3
        BootstrapExplainedAll(i,b)=(std(all_score(RandomSampleIndicesAll,i))^2)/(sum(std(data(RandomSampleIndicesAll,:)).^2));
        BootstrapExplainedBlind(i,b)=(std(all_score(RandomSampleIndicesBlind,i))^2)/(sum(std(data(RandomSampleIndicesBlind,:)).^2));
        BootstrapExplainedSight(i,b)=(std(all_score(RandomSampleIndicesSight,i))^2)/(sum(std(data(RandomSampleIndicesSight,:)).^2));
    end
    
end % across Bootstraps

stdBootstrapExplainedAll=std(BootstrapExplainedAll,0,2);
stdBootstrapExplainedBlind=std(BootstrapExplainedBlind,0,2);
stdBootstrapExplainedSight=std(BootstrapExplainedSight,0,2);


fprintf('Percent variance explained by each component of the final\n');
fprintf(' dimensions used (±SEM):\n\n');
fprintf('dim    sight    blind    both\n');
for i=1:3
    Outline=[num2str(i) '   ' ...
        num2str((std(all_score(indexsight,i))^2)/(sum(std(data(indexsight,:)).^2)),'%1.2f') ' ± ' ...
        num2str(stdBootstrapExplainedSight(i),'%1.2f') '  ' ...
        num2str((std(all_score(indexblind,i))^2)/(sum(std(data(indexblind,:)).^2)),'%1.2f') ' ± ' ...
        num2str(stdBootstrapExplainedBlind(i),'%1.2f') '  ' ...
        num2str((std(all_score(:,i))^2)/(sum(std(data(:,:)).^2)),'%1.2f') ' ± ' ...
        num2str(stdBootstrapExplainedAll(i),'%1.2f') ];
    
    fprintf([Outline '\n']);
end
Outline=['Total   ' ...
        num2str(sum(mean(BootstrapExplainedSight,2)),'%1.2f') ' ± ' ...
        num2str(std(sum(BootstrapExplainedSight,1)),'%1.2f') '  ' ...
        num2str(sum(mean(BootstrapExplainedBlind,2)),'%1.2f') ' ± ' ...
        num2str(std(sum(BootstrapExplainedBlind,1)),'%1.2f') '  ' ...
        num2str(sum(mean(BootstrapExplainedAll,2)),'%1.2f') ' ± ' ...
        num2str(std(sum(BootstrapExplainedAll,1)),'%1.2f') ];
 fprintf([Outline '\n']);
fprintf('\n\n');






% Report the Cohen's d effect size of blind vs sighted along each dimension

fprintf('Cohen-s d effect sizes for blind vs. sighted along each PC:\n');
for i=1:3
    d=(mean(all_score(indexsight,i))-mean(all_score(indexblind,i)))/ ( ( std(all_score(indexsight,i)) + std(all_score(indexblind,i)) ) /2);
    fprintf(['PC' num2str(i) ': ' num2str(d,'%1.1f') '\n']);
end
fprintf('\n\n');

% fprintf('Number of SEMs between blind vs. sighted along each PC:\n');
% for i=1:3
%     d=(mean(all_score(indexsight,i))-mean(all_score(indexblind,i)))/ ( ( std(all_score(indexsight,i))/sqrt(length(indexsight)) + std(all_score(indexblind,i))/sqrt(length(indexblind)) ) /2);
%     fprintf(['PC' num2str(i) ': ' num2str(d) '\n']);
% end
% fprintf('\n\n');





% Create horizontal bar plots of measurement weights for first three PCs
%  with ±1 SEM on the weights for the PCA on all subjects


YRanges=[-1 1; -1 1; -1 1];

figtmp=figure('name','Coefficients of PCs 1-3 derived from all subjects ±SEM');
titles={'PC1 - all subjects ±1SEM' 'PC2 - all subjects ±1SEM' 'PC3 - all subjects ±1SEM'};

for i=1:3
    subplot(3,1,i);
    hBars=bar(coeff(:,i),.05,'FaceColor',[.3 .2 .5],'LineStyle','none');
    hBaseline = get(hBars(1),'BaseLine');
    set(hBaseline,'LineStyle','--',...
        'Color',[0.3 0.3 0.3],...
        'LineWidth',0.5);
    hold on
    h=errorbar(coeff(:,i),sem_all_coeff(:,i),'Color',[0.8 0.8 0.8],'LineWidth',2);
    removeErrorBarEnds(h);
    scatter(linspace(1,NumMeasures,NumMeasures),coeff(:,i),40,[.5 0 .75],'+');
    hold on
    ylim(YRanges(i,:));
    pbaspect([2 2 1])
    title([titles{i}]);
    box off;
    hold off
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;








% Plot the sighted and blind subjects for PC1xPC2, color code the blind subgroups

figtmp=figure('name','All subjects PC1xPC2');
plot(all_score(indexsight,1),all_score(indexsight,2),'.k');
hold on
plot(all_score(indexanophthalmic,1),all_score(indexanophthalmic,2),'.r');
plot(all_score(indexcongenital,1),all_score(indexcongenital,2),'.g');
plot(all_score(indexpostnatal,1),all_score(indexpostnatal,2),'.b');
plot(all_score(indexlca,1),all_score(indexlca,2),'.y');
ylim([-4 4]);
xlim([-4 5]);
axis equal
hold off
box off;
title('All subjects PC1xPC2');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


% Now make a figure with just ellipses for the sub groups

figtmp=figure('name','All subjects PC1xPC2 - means ± SEM');

% Add a ±1 SD ellipse for all sighted and all blind

ellipse_size = 1; % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,2)],ellipse_size,'-k');
error_ellipse([all_score(indexblind,1),all_score(indexblind,2)],ellipse_size,'-r');

% Add ±1 SEM ellipses

ellipse_size = 1/sqrt(length(indexsight)); % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,2)],ellipse_size,'-k');

ellipse_size = 1/sqrt(length(indexblind)); % in units of SD
error_ellipse([all_score(indexblind,1),all_score(indexblind,2)],ellipse_size,'-r');

ellipse_size = 1/sqrt(length(indexanophthalmic)); % in units of SD
error_ellipse([all_score(indexanophthalmic,1),all_score(indexanophthalmic,2)],ellipse_size,'-r');

ellipse_size = 1/sqrt(length(indexcongenital)); % in units of SD
error_ellipse([all_score(indexcongenital,1),all_score(indexcongenital,2)],ellipse_size,'-g');

ellipse_size = 1/sqrt(length(indexpostnatal)); % in units of SD
error_ellipse([all_score(indexpostnatal,1),all_score(indexpostnatal,2)],ellipse_size,'-b');

ellipse_size = 1/sqrt(length(indexlca)); % in units of SD
error_ellipse([all_score(indexlca,1),all_score(indexlca,2)],ellipse_size,'-y');

ellipse_size = 1/sqrt(length(indexpostnatal_beforepuberty)); % in units of SD
error_ellipse([all_score(indexpostnatal_beforepuberty,1),all_score(indexpostnatal_beforepuberty,2)],ellipse_size,':r');

ellipse_size = 1/sqrt(length(indexpostnatal_afterpuberty)); % in units of SD
error_ellipse([all_score(indexpostnatal_afterpuberty,1),all_score(indexpostnatal_afterpuberty,2)],ellipse_size,':k');


plot(all_score(indexsight,1),all_score(indexsight,2),'.w');
plot(all_score(indexblind,1),all_score(indexblind,2),'.w');


ylim([-4 4])
xlim([-4 5])
axis equal
box off;
hold off
title('All subjects PC1xPC2 - means ± SEM');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;






% Plot the sighted and blind subjects for PC1xPC3, color code the blind subgroups

figtmp=figure('name','All subjects PC1xPC3');
plot(all_score(indexsight,1),all_score(indexsight,3),'.k');
hold on
plot(all_score(indexanophthalmic,1),all_score(indexanophthalmic,3),'.r');
plot(all_score(indexcongenital,1),all_score(indexcongenital,3),'.g');
plot(all_score(indexpostnatal,1),all_score(indexpostnatal,3),'.b');
plot(all_score(indexlca,1),all_score(indexlca,3),'.y');
ylim([-4 4]);
xlim([-4 4]);
axis equal
box off;
hold off
title('All subjects PC1xPC3');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


figtmp=figure('name','All subjects PC1xPC3 - means ± SEM');

% Add a ±1 SD ellipse for all sighted and all blind

ellipse_size = 1; % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,3)],ellipse_size,'-k');
error_ellipse([all_score(indexblind,1),all_score(indexblind,3)],ellipse_size,'-r');

% Add ±1 SEM ellipses

ellipse_size = 1/sqrt(length(indexsight)); % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,3)],ellipse_size,'-k');

ellipse_size = 1/sqrt(length(indexblind)); % in units of SD
error_ellipse([all_score(indexblind,1),all_score(indexblind,3)],ellipse_size,'-r');

ellipse_size = 1/sqrt(length(indexanophthalmic)); % in units of SD
error_ellipse([all_score(indexanophthalmic,1),all_score(indexanophthalmic,3)],ellipse_size,'-r');

ellipse_size = 1/sqrt(length(indexcongenital)); % in units of SD
error_ellipse([all_score(indexcongenital,1),all_score(indexcongenital,3)],ellipse_size,'-g');

ellipse_size = 1/sqrt(length(indexpostnatal)); % in units of SD
error_ellipse([all_score(indexpostnatal,1),all_score(indexpostnatal,3)],ellipse_size,'-b');

ellipse_size = 1/sqrt(length(indexlca)); % in units of SD
error_ellipse([all_score(indexlca,1),all_score(indexlca,3)],ellipse_size,'-y');

ellipse_size = 1/sqrt(length(indexpostnatal_beforepuberty)); % in units of SD
error_ellipse([all_score(indexpostnatal_beforepuberty,1),all_score(indexpostnatal_beforepuberty,2)],ellipse_size,':r');

ellipse_size = 1/sqrt(length(indexpostnatal_afterpuberty)); % in units of SD
error_ellipse([all_score(indexpostnatal_afterpuberty,1),all_score(indexpostnatal_afterpuberty,3)],ellipse_size,':k');

plot(all_score(indexsight,1),all_score(indexsight,3),'.w');
plot(all_score(indexblind,1),all_score(indexblind,3),'.w');

ylim([-4 4])
xlim([-4 4])
axis equal
box off;
hold off
title('All subjects PC1xPC3 - means ± SEM');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;



% Report the mean along each PC dimension for each group

fprintf('Group mean for sighted, anop, congen, postnata, lca\n');
for i=1:3
    Outline=['PC' num2str(i) ': '];
    Outline=[Outline num2str(mean(all_score(indexsight,i)),'%1.2g') ','];
    Outline=[Outline num2str(mean(all_score(indexanophthalmic,i)),'%1.2g') ','];
    Outline=[Outline num2str(mean(all_score(indexcongenital,i)),'%1.2g') ','];
    Outline=[Outline num2str(mean(all_score(indexpostnatal,i)),'%1.2g') ','];
    Outline=[Outline num2str(mean(all_score(indexlca,i)),'%1.2g') ','];
    fprintf([Outline '\n']);
end
fprintf('\n\n');



% Create bar plots of the group means and SEM for each dimension

figtmp=figure('name','Group Means ± SEM each dimension');
titles={'dim1' 'dim2' 'dim3'};


PCLabels={'anoph' 'congen' 'LCA' 'postnatal' 'sighted'};


for i=1:3
    subplot(3,1,i);
    tmpGroupDimScoreMeans=[mean(all_score(indexanophthalmic,i)) mean(all_score(indexcongenital,i)) mean(all_score(indexlca,i)) mean(all_score(indexpostnatal,i)) mean(all_score(indexsight,i))];
    tmpGroupDimScoreSEM=[std(all_score(indexanophthalmic,i))/sqrt(length(indexanophthalmic)) std(all_score(indexcongenital,i))/sqrt(length(indexcongenital)) std(all_score(indexlca,i))/sqrt(length(indexlca)) std(all_score(indexpostnatal,i))/sqrt(length(indexpostnatal)) std(all_score(indexsight,i))/sqrt(length(indexsight))];
    bar([1 2 3 4 5],tmpGroupDimScoreMeans);
    hold on
    errorbar([1 2 3 4 5],tmpGroupDimScoreMeans,tmpGroupDimScoreMeans,'o');
    hold on
    xlim([0 6]);
    ylim([-3 2]);
    set(gca,'Xtick',1:3,'XTickLabel',PCLabels);
    pbaspect([2 1 1])
    box off;
    title([titles{i}]);
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;
















% Test if there are differences in standard deviation in the PC1, PC3, and
% PC3 dimensions for the all subjects PCA.

fprintf('Test for differences in the standard deviation of the PC\n');
fprintf('  components between blind and sighted groups\n');

for i=1:3
    [p,stats]=vartestn([all_score(indexsight,i),all_score(indexblind,i)],'Display','off');
    fprintf(['PC',num2str(i),'- c2(1, N = ' num2str(length(indexsight)+length(indexblind)) ') = ',num2str(stats.chisqstat,'%1.2f'),', p=',num2str(p)],'%1.2g');
    fprintf('\n');
end
fprintf('\n\n');







% Remove the subjects identified as
% outliers in the initial box and whisker plots.

all_score_outlierclear=all_score;

for i=1:NumMeasures
    Outliers=find(abs(data(:,i))>OutThresh);
    if ~isempty(Outliers)
        all_score_outlierclear(Outliers,:)=NaN;
        fprintf(['Removed ' num2str(length(Outliers)) ' outliers\n\n']);
    end
end



% Report correlations of scores on the PC axes with age, gender, and
% supertentorial volume

fprintf('Correlation of scores on each PC axis with age, gender, and \n');
fprintf('  supertentorial volume:\n\n');

for i=1:3
    fprintf(['PC' num2str(i) ': ' num2str(nancorr(all_score_outlierclear(:,i),ages(:)),'%1.2f') '  ' num2str(nancorr(all_score_outlierclear(:,i),gender(:)),'%1.2f') '  ' num2str(nancorr(all_score_outlierclear(:,i),super_tentorial_volume(:)),'%1.2f') '\n']);
end

fprintf('\n\n');



%%%%%%%%%%%%%%%%%%%%%%%
% Time to work on the "other measures"
%%%%%%%%%%%%%%%%%%%%%%%







LabelsPC={'PC1' 'PC2' 'PC3'};
LabelsMeasures={'CBF' 'Sentences' 'FA'};





% Report the means, and t-test of the sighted vs. the blind

[~,p_vals,~,stats] = ttest2(othermeasures(indexsight,:),othermeasures(indexblind,:));
fprintf('Compare blind and all sighted on the CBF, Sentences, and FA measures.\n');
fprintf('Measure: [mean sighted , mean blind], t-test.\n');
for i=1:3
    Outline=[LabelsMeasures{i} ': [' num2str(nanmean(othermeasures(indexsight,i)),'%1.2f') ', ' num2str(nanmean(othermeasures(indexblind,i)),'%1.2f') ']'];
    Outline=[Outline ', t(' num2str(stats.df(i)) ') = ' num2str(stats.tstat(i),'%1.2f') ', p = ' num2str(p_vals(i),'%1.2g') '\n'];
    fprintf(Outline);
end
fprintf('\n\n');

[~,p_vals,~,stats] = ttest2(othermeasures(indexsight,:),othermeasures(indexcongenital,:));
fprintf('Compare early blind and all sighted on the CBF, Sentences, and FA measures.\n');
fprintf('Measure: [mean sighted , mean early blind], t-test.\n');
for i=1:3
    Outline=[LabelsMeasures{i} ': [' num2str(nanmean(othermeasures(indexsight,i)),'%1.2f') ', ' num2str(nanmean(othermeasures(indexcongenital,i)),'%1.2f') ']'];
    Outline=[Outline ', t(' num2str(stats.df(i)) ') = ' num2str(stats.tstat(i),'%1.2f') ', p = ' num2str(p_vals(i),'%1.2g') '\n'];
    fprintf(Outline);
end
fprintf('\n\n');

[~,p_vals,~,stats] = ttest2(othermeasures(indexsight,:),othermeasures(indexpostnatal,:));
fprintf('Compare late blind and all sighted on the CBF, Sentences, and FA measures.\n');
fprintf('Measure: [mean sighted , mean late blind], t-test.\n');
for i=1:3
    Outline=[LabelsMeasures{i} ': [' num2str(nanmean(othermeasures(indexsight,i)),'%1.2f') ', ' num2str(nanmean(othermeasures(indexpostnatal,i)),'%1.2f') ']'];
    Outline=[Outline ', t(' num2str(stats.df(i)) ') = ' num2str(stats.tstat(i),'%1.2f') ', p = ' num2str(p_vals(i),'%1.2g') '\n'];
    fprintf(Outline);
end
fprintf('\n\n');




% Report the correlation of each all-subject PC with the other measures,
% broken by sighted and blind subjects. Remove the subjects identified as
% outliers in the initial box and whisker plots.

all_score_outlierclear=all_score;
OutThresh=4.5;

for i=1:NumMeasures
    Outliers=find(abs(data(:,i))>OutThresh);
    if ~isempty(Outliers)
        all_score_outlierclear(Outliers,:)=NaN;
    end
end

LabelsPC={'PC1' 'PC2' 'PC3'};
LabelsMeasures={'CBF' 'Sentences' 'FA'};









% Conduct a multivariate regression for all subjects

d=3; % three data variables to model (CSF, sentences, FA)
p=3; % the three anatomical variation dimensions
n=NumSubjects; % NumSubjects rows, although most are NaNs

% Use a common ranges to plot these results

YRanges=[0 2.5; -1 3; 0.3 0.6];
XRanges=[-0.25 0.25; -.75 .75; -0.05 0.05];

% Create the matrix of data (Y) and matrix of predictors (X)

Y=[othermeasures(:,1) othermeasures(:,2) othermeasures(:,3)];
for i=1:d
    Y(:,i)=Y(:,i)-nanmean(Y(:,i));
end

X=[all_score_outlierclear(:,1) all_score_outlierclear(:,2) all_score_outlierclear(:,3)];
for i=1:p
    X(:,i)=X(:,i)-nanmean(X(:,i));
end

% Reformat X as a cell array to satisfy mvregress

CellX= repmat({zeros(d,p)},n,1);
for i=1:n
    for j = 1:d
        CellX{i} = kron(eye(d),X(i,:));
    end
end

% Regress

[beta,~,one_pop_model_E,V,one_pop_model_logL] = mvregress(CellX,Y);

B = reshape(beta,p,d)';
SEM = reshape(sqrt(diag(V)),p,d)';







% Report the beta weights on the principle components and SEM of these


fprintf('Beta weights (SEM) all\n')
fprintf('       PC1     PC2    PC3\n')
for i=1:3
    Outline=[LabelsMeasures{i} '   '];
    for k=1:3
        Outline=[Outline num2str(B(i,k),'%1.2g') '(' num2str(SEM(i,k),'%1.2g') ')  '];
    end
    fprintf(Outline);
    fprintf('\n');
end
fprintf('\n');

% Create a plot of this

figtmp=figure('name','Weights ±SEM on the components');
titles={'CBF' 'Sentences' 'FA'};
PCLabels={'PC1' 'PC2' 'PC3'};
PCPlotYRages=[-.25 .25; -.25 .25; -.025 .025];


for i=1:3
    subplot(3,1,i);
    bar([1 2 3],B(i,1:3));
    hold on
    errorbar([1 2 3],B(i,1:3),SEM(i,1:3),'o');
    hold on
    xlim([0 4]);
    ylim(PCPlotYRages(i,:));
    set(gca,'Xtick',1:3,'XTickLabel',PCLabels);
    pbaspect([2 1 1])
    box off;
    title([titles{i}]);
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;







% Create three scatter plots of model fit vs. measure

figtmp=figure('name','Correlation with other measures');

titles={'Fit vs CBF All' 'Fit vs Sentences All' 'Fit vs FA All'};

for i=1:3
    subplot(3,1,i);
    plot((B(i,1:3)*X(:,1:3)'),othermeasures(:,i),'.k');
    r=nancorr((B(i,1:3)*X(:,1:3)'),othermeasures(:,i));
    hold on
    lsline
    ylim(YRanges(i,:));
    xlim(XRanges(i,:));
    pbaspect([2 1 1])
    box off;
    title([titles{i} '  r=' num2str(r)]);
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


% Plot again with separate colors for sighted, congenital and postnatal

figtmp=figure('name','Correlation with other measures');

titles={'Fit vs CBF All' 'Fit vs Sentences All' 'Fit vs FA All'};

for i=1:3
    subplot(3,1,i);
    tmpfit=(B(i,1:3)*X(:,1:3)');
    plot(tmpfit(indexsight),othermeasures(indexsight,i),'.k');
    hold on
    plot(tmpfit(indexcongenital),othermeasures(indexcongenital,i),'.g');
    hold on
    plot(tmpfit(indexpostnatal),othermeasures(indexpostnatal,i),'.b');
    hold on
    lsline
    ylim(YRanges(i,:));
    xlim(XRanges(i,:));
    box off;
    pbaspect([2 1 1])
    title([titles{i}]);
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;





% Obtain the mean and SEM of the non-MPRAGE measures, broken down
%  by sub-group. Create a plot of this.

LabelsMeasures={'CBF' 'Sentences' 'FA'};
LabelsGroups={'Sight' 'Postnatal' 'Congen'};

OtherMeasureMean=zeros(3,3);
OtherMeasureSEM=zeros(3,3);


for i=1:3 % Loop across measures
    OtherMeasureMean(i,1)=nanmean(othermeasures(indexsight,i));
    OtherMeasureSEM(i,1)=nanstd(othermeasures(indexsight,i)) / sqrt(length(indexsight));
    
    OtherMeasureMean(i,2)=nanmean(othermeasures(indexpostnatal,i));
    OtherMeasureSEM(i,2)=nanstd(othermeasures(indexpostnatal,i)) / sqrt(length(indexpostnatal));
    
    OtherMeasureMean(i,3)=nanmean(othermeasures(indexcongenital,i));
    OtherMeasureSEM(i,3)=nanstd(othermeasures(indexcongenital,i)) / sqrt(length(indexcongenital));
end

figtmp=figure('name','Mean of other measures by group ±SEM');
titles={'CBF' 'Sentences' 'FA'};

for i=1:3
    subplot(3,1,i);
    errorbar([1 2 3],OtherMeasureMean(i,:),OtherMeasureSEM(i,:),'o');
    hold on
    ylim(YRanges(i,:));
    xlim([1 3]);
    set(gca,'Xtick',1:3,'XTickLabel',LabelsGroups);
    pbaspect([2 1 1])
    box off;
    title([titles{i}]);
end

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;













% Calculate F tests associated with the three PCs modeling the other
% measures

fprintf('F test of the 3 PC model of the measures for all subjects:\n');
for i=1:3
    [~,~,~,~,tmp_stats]=regress(Y(:,i),[X ones(1,NumSubjects)']);
    denom_f_df=NumSubjects-sum(isnan(Y(:,i)))-3;
    Outline=[titles{i} '- F (3,' num2str(denom_f_df) ') = ' num2str(tmp_stats(2),'%1.2f') ', p = ' num2str(tmp_stats(3),'%1.2g') '\n'];
    fprintf(Outline);
end
fprintf('\n\n');
















% Post-hoc: Make a plot of each group along a "VBI" type axis


[~,combined_data_score,~,~,~,~] = loo_pca(data);


figtmp=figure('name','Group Means ± SEM along a VBI');


PCLabels={'anoph' 'congen' 'LCA' 'postnatal' 'sighted'};


    tmpGroupDimScoreMeans=[mean(combined_data_score(indexanophthalmic,1)) mean(combined_data_score(indexcongenital,1)) mean(combined_data_score(indexlca,1)) mean(combined_data_score(indexpostnatal,1)) mean(combined_data_score(indexsight,1))];
    tmpGroupDimScoreSEM=[std(combined_data_score(indexanophthalmic,1))/sqrt(length(indexanophthalmic)) std(combined_data_score(indexcongenital,1))/sqrt(length(indexcongenital)) std(combined_data_score(indexlca,1))/sqrt(length(indexlca)) std(combined_data_score(indexpostnatal,1))/sqrt(length(indexpostnatal)) std(combined_data_score(indexsight,1))/sqrt(length(indexsight))];
    bar([1 2 3 4 5],tmpGroupDimScoreMeans);
    hold on
    errorbar([1 2 3 4 5],tmpGroupDimScoreMeans,tmpGroupDimScoreSEM,'o');
    hold on
    xlim([0 6]);
    ylim([-3 2]);
    set(gca,'Xtick',1:3,'XTickLabel',PCLabels);
    pbaspect([2 1 1])
    box off;

saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;











% Post-hoc: Check if there is a difference in correlation of the early
% blind compared to the late blind / sighted of anatomical dimension 3 and
% the optic radiation / splenium FA

fprintf(['Corr of early blind FA with anatomy: ' num2str(nancorr((B(3,1:3)*X(indexcongenital,1:3)'),othermeasures(indexcongenital,3)),'%1.2f') '\n' ]);
fprintf(['Corr of late blind and sighted FA with anatomy: ' num2str(nancorr((B(3,1:3)*X([indexpostnatal indexsight],1:3)'),othermeasures([indexpostnatal indexsight],3)),'%1.2f') '\n' ]);
fprintf(['For early corr, n=' num2str(length(othermeasures([indexsight indexpostnatal],3))-sum(isnan((othermeasures([indexsight indexpostnatal],3))))) '\n']);
fprintf(['For late and sighted corr, n=' num2str(length(othermeasures(indexcongenital,3))-sum(isnan((othermeasures(indexcongenital,3))))) '\n']);
fprintf('\n\n');
fprintf(['Corr of blind BOLD fMRI with anatomy: ' num2str(nancorr((B(2,1:3)*X([indexpostnatal indexcongenital],1:3)'),othermeasures([indexpostnatal indexcongenital],2)),'%1.2f') '\n' ]);
fprintf(['Corr of sighted BOLD fMRI with anatomy: ' num2str(nancorr((B(2,1:3)*X(indexsight,1:3)'),othermeasures(indexsight,2)),'%1.2f') '\n' ]);
fprintf(['For early and late corr, n=' num2str(length(othermeasures([indexcongenital indexpostnatal],3))-sum(isnan((othermeasures([indexcongenital indexpostnatal],3))))) '\n']);
fprintf(['For sighted corr, n=' num2str(length(othermeasures(indexsight,3))-sum(isnan((othermeasures(indexsight,3))))) '\n']);
fprintf('\n\n');




% Post-hoc: Give the p-value for categorical covariates that model
%  blind v sighted, congenital vs. postnatal, and [anop vs. congent + LCA
%  vs. postnatal] for each dimension

x1=zeros(118,1);
x1(indexsight)=1;
x1(indexblind)=-1;

x2=zeros(118,1);
x2([indexlca indexpostnatal])=1;
x2([indexanophthalmic indexcongenital])=-1;

x3=zeros(118,1);
x3([indexpostnatal indexcongenital])=1;
x3([indexlca indexanophthalmic])=-1;

X=[x1 x2 x3];

for i=1:3
    [betas,dev,stats] = glmfit(X,all_score(:,i));
    stats.t
    stats.p
end


% Post-hoc: See if there is a correlation in the post-natally blind of Dim1
% and years of blindness:

fprintf(['Correlation of years of max blindness with Dim1 score, postnatal group: ' num2str(nancorr(ages(indexpostnatal)-blindness_onset_max_years(indexpostnatal,2),all_score_outlierclear(indexpostnatal,1)),'%1.2f') '\n\n']);


% Post-hoc: See if there is a difference in mean anatomical scores on each
% dimension between the sighted and those with blindness onset after the
% age of 25

[h,p,ci,stats]=ttest2(all_score(indexsight,1),all_score(indexpostnatal_afterpuberty,1));
fprintf([ 't (' num2str(length(indexsight)+length(indexpostnatal_afterpuberty)-2) ')=' num2str(stats.tstat,'%1.2f') ', p=' num2str(p,'%1.2g') '\n' ]);
[h,p,ci,stats]=ttest2(all_score(indexsight,2),all_score(indexpostnatal_afterpuberty,2));
fprintf([ 't (' num2str(length(indexsight)+length(indexpostnatal_afterpuberty)-2) ')=' num2str(stats.tstat,'%1.2f') ', p=' num2str(p,'%1.2g') '\n' ]);
[h,p,ci,stats]=ttest2(all_score(indexsight,3),all_score(indexpostnatal_afterpuberty,3));
fprintf([ 't (' num2str(length(indexsight)+length(indexpostnatal_afterpuberty)-2) ')=' num2str(stats.tstat,'%1.2f') ', p=' num2str(p,'%1.2g') '\n' ]);

[h,p,ci,stats]=ttest2(all_score(indexpostnatal_beforepuberty,1),all_score(indexpostnatal_afterpuberty,1));
fprintf([ 't (' num2str(length(indexpostnatal_beforepuberty)+length(indexpostnatal_afterpuberty)-2) ')=' num2str(stats.tstat,'%1.2f') ', p=' num2str(p,'%1.2g') '\n' ]);
[h,p,ci,stats]=ttest2(all_score(indexpostnatal_beforepuberty,2),all_score(indexpostnatal_afterpuberty,2));
fprintf([ 't (' num2str(length(indexpostnatal_beforepuberty)+length(indexpostnatal_afterpuberty)-2) ')=' num2str(stats.tstat,'%1.2f') ', p=' num2str(p,'%1.2g') '\n' ]);
[h,p,ci,stats]=ttest2(all_score(indexpostnatal_beforepuberty,3),all_score(indexpostnatal_afterpuberty,3));
fprintf([ 't (' num2str(length(indexpostnatal_beforepuberty)+length(indexpostnatal_afterpuberty)-2) ')=' num2str(stats.tstat,'%1.2f') ', p=' num2str(p,'%1.2g') '\n' ]);


% Post-hoc: Make a figure that has the coefficients for the anatomical
% dimensions defined using the all group data with no rotation.

fig_s4_coeff=pca(data);
[~,fig_s4_score,~,~,~,~] = loo_pca(data);
[ boot_fig_s4_coeff, sem_fig_s4_coeff, ~, ~ ] = boot_pca(data);



% Create a matrix plot of the loading on the significant anatomical
%  dimensions in the first three principal components

M=fig_s4_coeff;
noonsig=find(abs(fig_s4_coeff./sem_fig_s4_coeff)<1.95);
M(noonsig)=0;

figtmp=figure('name','Anatomical weights matrix, significant values only, coeffs for combo data');
h = imagesc(M);
axis square;
caxis([-0.7 0.7]);
C=polarmap();
colormap(C);
colorbar;
set(gca, 'visible', 'off') ;
title('Anatomical weights on 3 dimensions - Fig S4');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;


% Create a matrix plot of the loading on all anatomical
%  dimensions in the first three principal components

M=fig_s4_coeff;

figtmp=figure('name','Anatomical weights matrix all values, coeffs for combo data');
h = imagesc(M);
axis square;
caxis([-0.7 0.7]);
C=polarmap();
colormap(C);
colorbar;
set(gca, 'visible', 'off') ;
title('Anatomical weights on 3 dimensions - Fig S4');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.png'], 'png');
FigIndex=FigIndex+1;



% Dump out the columns for Table S4

for j=1:9
Outline=' ';
for i=1:3
Outline=[Outline num2str(coeff(j,i),'%1.2g') '±' num2str(sem_all_coeff(j,i),'%1.2g') '  '];
end
fprintf([Outline '\n']);
end



% Drop a footer to mark out the end of BlindnessNOS text output

fprintf('\n\n\n***************************************************\n\n')


fprintf('\n\nLCA1 ANALYSIS\n\n');


% Print a table of mean and SEM of the anatomical measures for each group

fprintf('Mean ± SEM for each measure for each group:\n');
fprintf('        Sighted     AllBlind  Anop        Congenital  Postnatal   AllLCA     LCA1        RPE65       CRB1        CEP290:\n');
 for i=1:NumMeasures
     Outline=['Measure ' num2str(i) '  ' ];
     Outline=[Outline num2str(mean(data(indexsight,i)),'%1.2f') '±' num2str(std(data(indexsight,i))/sqrt(length(indexsight)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(data(indexblind,i)),'%1.2f') '±' num2str(std(data(indexblind,i))/sqrt(length(indexblind)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(data(indexanophthalmic,i)),'%1.2f') '±' num2str(std(data(indexanophthalmic,i))/sqrt(length(indexanophthalmic)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(data(indexcongenital,i)),'%1.2f') '±' num2str(std(data(indexcongenital,i))/sqrt(length(indexcongenital)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(data(indexpostnatal,i)),'%1.2f') '±' num2str(std(data(indexpostnatal,i))/sqrt(length(indexpostnatal)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(data(indexlca,i)),'%1.2f') '±' num2str(std(data(indexlca,i))/sqrt(length(indexlca)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(data(indexlca1,i)),'%1.2f') '±' num2str(std(data(indexlca1,i))/sqrt(length(indexlca1)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(data(indexrpe65,i)),'%1.2f') '±' num2str(std(data(indexrpe65,i))/sqrt(length(indexrpe65)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(data(indexcrb1,i)),'%1.2f') '±' num2str(std(data(indexcrb1,i))/sqrt(length(indexcrb1)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(data(indexcep290,i)),'%1.2f') '±' num2str(std(data(indexcep290,i))/sqrt(length(indexcep290)),'%1.2f') '  '];
     fprintf([Outline '\n']);
 end
 fprintf('\n\n');
 
 
 
% Plot the sighted and blind subjects and LCA subjects for PC1xPC2, color code the blind subgroups

figtmp=figure('name','LCA subjects PC1xPC2');
plot(all_score(indexlca1,1),all_score(indexlca1,2),'.m');
hold on
plot(all_score(indexrpe65,1),all_score(indexrpe65,2),'.g');
plot(all_score(indexcrb1,1),all_score(indexcrb1,2),'.b');
plot(all_score(indexcep290,1),all_score(indexcep290,2),'.y');
ylim([-4 4]);
xlim([-4 5]);
axis equal
hold off
box off;
title('All subjects PC1xPC2');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


% Now make a figure with just ellipses for the sub groups

figtmp=figure('name','All subjects PC1xPC2 - means ± SEM');

% Add a ±1 SD ellipse for all sighted and all blind

ellipse_size = 1; % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,2)],ellipse_size,'-k');
error_ellipse([all_score(indexblind,1),all_score(indexblind,2)],ellipse_size,'-r');

% Add ±1 SEM ellipses

ellipse_size = 1/sqrt(length(indexsight)); % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,2)],ellipse_size,'-k');

ellipse_size = 1/sqrt(length(indexblind)); % in units of SD
error_ellipse([all_score(indexblind,1),all_score(indexblind,2)],ellipse_size,'-r');

ellipse_size = 1/sqrt(length(indexlca1)); % in units of SD
error_ellipse([all_score(indexlca1,1),all_score(indexlca1,2)],ellipse_size,'-m');

ellipse_size = 1/sqrt(length(indexrpe65)); % in units of SD
error_ellipse([all_score(indexrpe65,1),all_score(indexrpe65,2)],ellipse_size,'-g');

ellipse_size = 1/sqrt(length(indexcrb1)); % in units of SD
error_ellipse([all_score(indexcrb1,1),all_score(indexcrb1,2)],ellipse_size,'-b');

ellipse_size = 1/sqrt(length(indexcep290)); % in units of SD
error_ellipse([all_score(indexcep290,1),all_score(indexcep290,2)],ellipse_size,'-y');


plot(all_score(indexsight,1),all_score(indexsight,2),'.w');
plot(all_score(indexblind,1),all_score(indexblind,2),'.w');


ylim([-4 4])
xlim([-4 5])
axis equal
box off;
hold off
title('All subjects PC1xPC2 - means ± SEM');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;



 
% Plot the sighted and blind subjects and LCA subjects for PC1xPC3, color code the blind subgroups

figtmp=figure('name','LCA subjects PC1xPC2');
plot(all_score(indexlca1,1),all_score(indexlca1,3),'.m');
hold on
plot(all_score(indexrpe65,1),all_score(indexrpe65,3),'.g');
plot(all_score(indexcrb1,1),all_score(indexcrb1,3),'.b');
plot(all_score(indexcep290,1),all_score(indexcep290,3),'.y');
ylim([-4 4]);
xlim([-4 4]);
axis equal
hold off
box off;
title('All subjects PC1xPC3');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


% Now make a figure with just ellipses for the sub groups

figtmp=figure('name','All subjects PC1xPC3 - means ± SEM');

% Add a ±1 SD ellipse for all sighted and all blind

ellipse_size = 1; % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,3)],ellipse_size,'-k');
error_ellipse([all_score(indexblind,1),all_score(indexblind,3)],ellipse_size,'-r');

% Add ±1 SEM ellipses

ellipse_size = 1/sqrt(length(indexsight)); % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,3)],ellipse_size,'-k');

ellipse_size = 1/sqrt(length(indexblind)); % in units of SD
error_ellipse([all_score(indexblind,1),all_score(indexblind,3)],ellipse_size,'-r');

ellipse_size = 1/sqrt(length(indexlca1)); % in units of SD
error_ellipse([all_score(indexlca1,1),all_score(indexlca1,3)],ellipse_size,'-m');

ellipse_size = 1/sqrt(length(indexrpe65)); % in units of SD
error_ellipse([all_score(indexrpe65,1),all_score(indexrpe65,3)],ellipse_size,'-g');

ellipse_size = 1/sqrt(length(indexcrb1)); % in units of SD
error_ellipse([all_score(indexcrb1,1),all_score(indexcrb1,3)],ellipse_size,'-b');

ellipse_size = 1/sqrt(length(indexcep290)); % in units of SD
error_ellipse([all_score(indexcep290,1),all_score(indexcep290,3)],ellipse_size,'-y');


plot(all_score(indexsight,1),all_score(indexsight,2),'.w');
plot(all_score(indexblind,1),all_score(indexblind,2),'.w');


ylim([-4 4])
xlim([-4 4])
axis equal
box off;
hold off
title('All subjects PC1xPC3 - means ± SEM');
saveas(figtmp, [OutFileStem 'BlindNOS_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;













% Report the mean along each PC dimension for each group

fprintf('Mean ± SEM for each dimension for each group:\n');
fprintf('        Sighted     AllBlind  Anop        Congenital  Postnatal   AllLCA     LCA1        RPE65       CRB1        CEP290:\n');
 for i=1:3
     Outline=['Measure ' num2str(i) '  ' ];
     Outline=[Outline num2str(mean(all_score(indexsight,i)),'%1.2f') '±' num2str(std(all_score(indexsight,i))/sqrt(length(indexsight)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexblind,i)),'%1.2f') '±' num2str(std(all_score(indexblind,i))/sqrt(length(indexblind)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexanophthalmic,i)),'%1.2f') '±' num2str(std(all_score(indexanophthalmic,i))/sqrt(length(indexanophthalmic)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexcongenital,i)),'%1.2f') '±' num2str(std(all_score(indexcongenital,i))/sqrt(length(indexcongenital)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexpostnatal,i)),'%1.2f') '±' num2str(std(all_score(indexpostnatal,i))/sqrt(length(indexpostnatal)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexlca,i)),'%1.2f') '±' num2str(std(all_score(indexlca,i))/sqrt(length(indexlca)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexlca1,i)),'%1.2f') '±' num2str(std(all_score(indexlca1,i))/sqrt(length(indexlca1)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexrpe65,i)),'%1.2f') '±' num2str(std(all_score(indexrpe65,i))/sqrt(length(indexrpe65)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexcrb1,i)),'%1.2f') '±' num2str(std(all_score(indexcrb1,i))/sqrt(length(indexcrb1)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexcep290,i)),'%1.2f') '±' num2str(std(all_score(indexcep290,i))/sqrt(length(indexcep290)),'%1.2f') '  '];
     fprintf([Outline '\n']);
 end
 fprintf('\n\n');

 


% Print a table of mean and SEM of the other measures for each group

fprintf('Mean ± SEM for each dimension for each group:\n');
fprintf('        Sighted     AllBlind  Anop        Congenital  Postnatal   AllLCA     LCA1        RPE65       CRB1        CEP290:\n');
 for i=1:3
     Outline=['Measure ' num2str(i) '  ' ];
     Outline=[Outline num2str(mean(all_score(indexsight,i)),'%1.2f') '±' num2str(std(all_score(indexsight,i))/sqrt(length(indexsight)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexblind,i)),'%1.2f') '±' num2str(std(all_score(indexblind,i))/sqrt(length(indexblind)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexanophthalmic,i)),'%1.2f') '±' num2str(std(all_score(indexanophthalmic,i))/sqrt(length(indexanophthalmic)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexcongenital,i)),'%1.2f') '±' num2str(std(all_score(indexcongenital,i))/sqrt(length(indexcongenital)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexpostnatal,i)),'%1.2f') '±' num2str(std(all_score(indexpostnatal,i))/sqrt(length(indexpostnatal)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexlca,i)),'%1.2f') '±' num2str(std(all_score(indexlca,i))/sqrt(length(indexlca)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexlca1,i)),'%1.2f') '±' num2str(std(all_score(indexlca1,i))/sqrt(length(indexlca1)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexrpe65,i)),'%1.2f') '±' num2str(std(all_score(indexrpe65,i))/sqrt(length(indexrpe65)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexcrb1,i)),'%1.2f') '±' num2str(std(all_score(indexcrb1,i))/sqrt(length(indexcrb1)),'%1.2f') '  '];
     Outline=[Outline num2str(mean(all_score(indexcep290,i)),'%1.2f') '±' num2str(std(all_score(indexcep290,i))/sqrt(length(indexcep290)),'%1.2f') '  '];
     fprintf([Outline '\n']);
 end
 fprintf('\n\n');

 

 
 
 
 % Drop a footer to mark out the end of LCA1 text output

fprintf('\n\n\n***************************************************\n\n')

close all;
