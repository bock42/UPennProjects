% Figures and Analyses for LCA1 paper
%  Much of the analysis begins with that used for the BlindnessNOS paper

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

% Hard coded ID codes for the LCA1 subjects

lca1_id=['A012813W';'A013013U';'B062813O';'C022813F';'G022513T';'G062713O'];


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
 
 


% Make a figure with  ellipses for the sub groups

figtmp=figure('name','All subjects PC1xPC2 - means ± SEM');

% Add a ±1 SD ellipse for all sighted and all blind

ellipse_size = 1; % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,2)],ellipse_size,'-k');
error_ellipse([all_score(indexblind,1),all_score(indexblind,2)],ellipse_size,'-r');

hold on

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

% Plot the individual LCA1 subejcts

plot(all_score(indexlca1,1),all_score(indexlca1,2),'.m');
text(all_score(indexlca1,1),all_score(indexlca1,2), lca1_id, 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right')


ylim([-4 4])
xlim([-4 5])
axis equal
box off;
hold off
title('All subjects PC1xPC2 - means ± SEM');
saveas(figtmp, [OutFileStem 'LCA1_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;




% Make a figure with  ellipses for the sub groups

figtmp=figure('name','All subjects PC1xPC3 - means ± SEM');

% Add a ±1 SD ellipse for all sighted and all blind

ellipse_size = 1; % in units of SD
error_ellipse([all_score(indexsight,1),all_score(indexsight,3)],ellipse_size,'-k');
error_ellipse([all_score(indexblind,1),all_score(indexblind,3)],ellipse_size,'-r');

hold on;

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

% Plot the individual LCA1 subejcts

plot(all_score(indexlca1,1),all_score(indexlca1,3),'.m');

plot(all_score(indexlca1,1),all_score(indexlca1,3),'.m');
text(all_score(indexlca1,1),all_score(indexlca1,3), lca1_id, 'VerticalAlignment','bottom', ...
                             'HorizontalAlignment','right')


ylim([-4 4])
xlim([-4 4])
axis equal
box off;
hold off
title('All subjects PC1xPC3 - means ± SEM');
saveas(figtmp, [OutFileStem 'LCA1_Fig' num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;










% Report the mean along each PC dimension for each group

fprintf('Mean ± SEM for each dimension for each group:\n');
fprintf('        Sighted     AllBlind  Anop        Congenital  Postnatal   AllLCA     LCA1        RPE65       CRB1        CEP290:\n');
 for i=1:3
     Outline=['Dim ' num2str(i) '  ' ];
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

%close all;
