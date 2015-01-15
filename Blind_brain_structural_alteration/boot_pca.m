function [ mean_coeff, sem_coeff, mean_explained, sem_explained ] = boot_pca( X, T, numdims )

% Implements a bootstrap resample to estimate error on the coefficients

RotateFlag=true;

if (~exist('T', 'var'))
    RotateFlag=false;
end

sizer=size(X);
IndicesToSample=linspace(1,sizer(1),sizer(1));


nBootstraps = 1000;

lengthIndicesToSample = sizer(1);

BootstrapCoeff = zeros(sizer(2),sizer(2),nBootstraps);
BootstrapExplained = zeros(sizer(2),nBootstraps);

for b = 1:nBootstraps
    
    
    % This will sample with replacement from the rows of X,
    % matching the total number of available rows in X.
    
    SampleWithReplace= ...
        datasample(linspace(1,lengthIndicesToSample, ...
        lengthIndicesToSample),lengthIndicesToSample,'Replace',true);
    
    RandomSampleIndices = IndicesToSample(SampleWithReplace);
    
    % Build a resampled X
    
    for i=1:lengthIndicesToSample
        if (i==1)
            ResampleX=[X(RandomSampleIndices(i),:)];
        else
            ResampleX=[ResampleX;X(RandomSampleIndices(i),:)];
        end % if the first row
    end % for number of samples
    
    [ TmpCoeff,~,~,~,explained,~ ] = pca( ResampleX );
    
    if RotateFlag
        TmpCoeff(:,1:numdims) = TmpCoeff(:,1:numdims)*T;
        TmpCoeff(:,numdims+1:end)= NaN;
    end
    
    BootstrapCoeff(:,:,b)=TmpCoeff;
    BootstrapExplained(:,b) = explained;
    
end % across Bootstraps



% Note that when we take the standard deviation across the bootstraps
% we are obtaining an estimate of the standard error of the mean in the
% full length data.

mean_coeff = mean(BootstrapCoeff,3);
sem_coeff = std(BootstrapCoeff,0,3);

mean_explained = mean(BootstrapExplained,2);
sem_explained = std(BootstrapExplained,0,2);

end

