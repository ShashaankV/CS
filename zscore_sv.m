%zscore with missing values, default over rows
function A=zscore_sv(G,miss_value,sub_tag,varargin)
    [m,n]=size(G);
    A=zeros(m,n);
    for i=1:n
        ind=G(:,i)~=miss_value;
        indnull=logical(1-ind);
        if sum(sum(G(:,1)-floor(G(:,1))))==0
            p=calc_maf(G(:,i),miss_value); %calculate p for current measure only
        end
        %calculate mean and standard deviation
        if nargin>3
                if strcmp(vargin{1},'T')
                    if sum(sum(G(:,1)-floor(G(:,1))))==0 %test if elements in the first column of G are integers, if true assume Binomial matrix, if not assume Gaussian
                        mu=2*p;
                        sig=sqrt(2*p*(1-p));
                    else
                        mu=mean(G(ind,i)); %for non-Binomial, calculate the empirical based on non-missing
                        sig=std(G(ind,i));
                    end
                end
        else
            mu=mean(G(ind,i)); %for non-Binomial, calculate the empirical based on non-missing
            sig=std(G(ind,i));
        end
        %zero-mean and scale variance with missing substitute depending on sub_tag
        if strcmp(sub_tag,'random') %substitute random value prior to normalization, then normalize all
            if sum(sum(G(:,1)-floor(G(:,1))))==0 %test if elements in the first column of G are integers, if true assume Binomial matrix, if not assume Gaussian
                x=binornd(2,p,[m,1]);
            else %if not Binomial assume Gaussian
                x=randn(mu,sig,[m,1]);
            end
            G(indnull,i)=x(indnull);
            if sig~=0
                A(ind,i)=(G(:,i)-mu)/sig;
            else
                A(ind,i)=0; %if variance is zero, then zero everything (due to order of limits), this may occur if the sample has no variance due to finite size, etc.
            end
        elseif strcmp(sub_tag,'zero') %normalize only non-missing then swap in zero for missing
            if sig~=0
                A(ind,i)=(G(ind,i)-mu)/sig;
            else
                A(ind,i)=0; %if variance is zero, then zero everything (due to order of limits)
            end
            A(indnull,i)=0; %zero missing
        end
    end
end