function [ coeff,score,latent,tsquared,explained,mu ] = loo_pca( X, T, numdims )

RotateFlag=true;

if (~exist('T', 'var'))
    RotateFlag=false;
end

% Implements a LOO approach to calculation of scores for the rows of X

sizer=size(X);

[ coeff,score,latent,tsquared,explained,mu ] = pca( X );
if RotateFlag
    coeff(:,1:numdims) = coeff(:,1:numdims)*T;
    coeff(:,numdims+1:end)= NaN;
end

loo_score=score.*0;

for i=1:sizer(1)
    indexer=linspace(1,sizer(1),sizer(1));
    indexer(i)=NaN;
    indexer=find(indexer > 0);
    subX=X(indexer,:);
    sub_coeff = pca( X );
    if RotateFlag
        sub_coeff(:,1:numdims) = sub_coeff(:,1:3)*T;
        sub_coeff(:,numdims+1:end)= NaN;
    end
    for component=1:sizer(2)
        loo_score(i,component)=sum(sub_coeff(:,component)'.*X(i,:));
    end
end


score=loo_score;


end

