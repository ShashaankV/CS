%References:
%1. Friedman,  J. H., Hastie, T., Hoefling H., and Tibshirani, R. Pathwise Coordinate Optimization. Annals Applied Statistics 1, 302-332 (2007)
%2. Friedman,  J. H., Hastie, T. and Tibshirani, R. Regularized Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1) (2010)

function varargout=lasso_sv_v2(G,miss_value,y,h2,varargin)
    [mg,ng]=size(G); %A is expected to be a matrix of standard normal columns, with subjects as rows (mg) and genetic variants as columns (ng)
    A=zscore_sv(G,miss_value,'zero')/sqrt(mg);
    keep=zeros(size(A,2));
    y=zscore(y)/sqrt(mg);
    if nargin>4
        if nargin==5
            fprintf('Not enough FE arguments')
            return
        else
            bool_fe=1; 
%             FE=zscore(varargin{1},[],1)/sqrt(mg);
            FE=zscore_sv(varargin{1},miss_value,'zero')/sqrt(mg); %comment out if using testing FE against known genotype vector (otherwise zscore w/o accounting for missing will deflate the fit) 
            A=cat(2,A,FE);
            FEflag=varargin{2};
            keep=zeros(size(A,2));
            if strcmp(FEflag,'keep')
                keep(ng+1:end)=1; %always keep active and never softhreshold
            end
        end
    else
        bool_fe=0;
    end
    
    lambda_max=max(abs(A'*y));
    %calc lambda min
        %use the median maximum value of A transpose times a noise vector z;
        %median taken from 1,000 samples of z (1,000 was arbirtrarily chosen)
        %z is scaled by sige (and 1/sqrt(ng) as above for A and y); also included is a
        %sampling error term 1/sqrt(ng)
    sige=sqrt(min(1,1-h2+1/mg)); %assumes that the variance of y=1
    z=sige*randn(mg,1000)/sqrt(mg);
    lambda_min=median(max(abs(A'*z)));
    %calc lambda series
    loghi = log(lambda_max);
    loglo = log(lambda_min);
    logrange = loghi - loglo;
    nLambda=100;
    interval = -logrange/(nLambda-1);
    lambda_ = exp(loghi:interval:loglo)';
    for lambi=1:length(lambda_)
        if lambi==1
            xhat=A'*y; %initialize coefficients by LSE if first lambda; otherwise, use fit from prior lambda
        end
        lambda=lambda_(lambi);
        r=y-A*xhat;
        kill=0; %kill variable for lasso subloop
        iter=0; %iteration count lasso subloop, used for indexing error
        err=[];
        %initialize active set, all
        active=ones(length(xhat),1); %active keeps track of current nonzeros; for iter==1, scan all parameters
        while kill==0
            iter=iter+1;
            ind=find(active~=0); %grab index of current nonzeros, and fit only these
            for ji=1:length(ind)
                j=ind(ji);
                xjold = xhat(j);
                xhat(j) = xjold+A(:,j)'*(r);
                if keep(j)==0
                    xhat(j) = sign(xhat(j)).*max(abs(xhat(j))-lambda,0); %if not in keep then softhreshold
                    active(j)=xhat(j)~=0; %if not in keep, then allow for deletion from active set 
                end
                r=r-A(:,j)*(xhat(j)-xjold); %update total residual by change in partial residual
            end
            ind=find(xhat~=0);
            err(iter)=sum(r.^2)+lambda*sum(abs(xhat(ind)));
            if iter>1
                if 1-min([err(iter-1),err(iter)])/max([err(iter-1),err(iter)])<10^-4  | isfinite(err(iter))~=1 %end lamba subloop if delta error is below threshold or if err==NaN
                    kill=1;
                end
            end
        end
    end
    if bool_fe==1
        if nargout==2
            varargout={xhat(1:ng),xhat(ng+1:end)};
        elseif nargout==1
            varargout={xhat(1:ng)};
        end
    else 
        varargout={xhat(1:ng)};
    end
end