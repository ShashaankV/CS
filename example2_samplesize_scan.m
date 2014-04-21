%Example script for running sample size scan
%correlated Normal random variables with chromosome 22 correlation structure
%Prior to running any examples, download the GenomicData contents to a separate (sub)folder and
%add this to the MATLAB path.

clc
clear

GDpath='./GenomicData'; %define the path for the GenomicData folder, the correlation matrix will be unpacked here
unpack_chr22_SNPcorrelations(GDpath)
clear

load('chr22_SNPcorrelationmatrix','R') %load correlation matrix or define it here 

n=length(R); %number of parameter 
m=floor(1.5*n); %some number of samples
miss=.05; %probability of a missing genotype
miss_value=-1; %indicator value for missing genotypes

G=G_Normal(m,n,miss,miss_value,R); %R is the desired correlation matrix for G

%define additive and coefficient model parameters
xs=125; %number of nonzeros
xsign='posneg'; %choices= pos, neg, posneg
xmaf='random'; %(only required for the Binomial genotype matrix; ignored otherwise) choices= maf_low, maf_high, random
xtype='Uniform'; %choices= Uniform, Hyperexponential1, Hyperexponential2
h2=.5;

[x,y]=y_additive(G,miss_value,h2,xs,xtype,xsign,xmaf); %the script accommodates Binomial and Normal matrices 

save('syn_data_ex2','G','miss_value','x','y','h2') 
clear
load('syn_data_ex2') %uncomment this line and comment the above, if reusing the synthetic data
                    %alternatively load your own data
                        %you will need to define G, y, miss_value, and h2...x is optional..if x is absent, then measures below need to be edited 

[mg,ng]=size(G);
mmax=8000; %define the maximum number of samples, subjects to scan
dm=floor(mmax/20); %survey 20 evenly spaced sample sizes
m_=dm:dm:mmax; %substitute your sampling scheme
%uncomment below if testing fixed effects using some of the known nonzeros
    %if using x nonzeros for FE, then lasso_sv should be edited to account for miss_value as noted in the function
% ind=find(x~=0);
% FE=G(:,ind(1:5));
% xFE=x(ind(1:5));
% G(:,ind(1:5))=[];
% x(ind(1:5))=[];
% FEflag='keep';
    
ind=0;
for m=m_
    ind=ind+1
    indm=randperm(mg,m); %sample random subset of subjects
    %uncomment next two lines if above fails due to older version of MATLAB 
%     indm=randperm(mg); 
%     indm=indm(1:m);
    xhat=lasso_sv(G(indm,:),miss_value,y(indm,1),h2);  
%     [xhat,xhatFE]=lasso_sv_fe(G(indm,:),miss_value,y(indm,1),h2,FE(indm,:),FEflag); %FE and FEflag are optional but need to be paired
                                                                                        %FE is assumed to have no missing values
    
    %next calculate a number of selection measures, note all but muP
    %require knowledge of the true positives
    %a set of proxies may be used instead of x for FPR and PPV by substituting a boolean vector 
    %where 1=acceptable nonzero (proxy) and 0=zero coefficient
    dat(ind,:)=[NE(x,xhat), FPR(x~=0,xhat~=0), PPV(x~=0,xhat~=0), medianP(xhat,G(indm,:),miss_value,y(indm,1))];
end

%plot results
fs=14;
figure(1)
clf
subplot(1,2,1)
[ax,h]=plotyy(m_,dat(:,1),m_,dat(:,3));
ylim(ax(1),[0 max(dat(:,1))])
ylim(ax(2),[0 1])
xlabel(ax(1),'sample size','fontsize',fs)
ylabel(ax(1),'NE','fontsize',fs)
ylabel(ax(2),'PPV','fontsize',fs)

subplot(1,2,2)
[ax,h]=plotyy(m_,dat(:,2),m_,dat(:,4));
xlabel(ax(1),'sample size','fontsize',fs)
ylabel(ax(1),'FPR','fontsize',fs)
ylabel(ax(2),'\mu_{P-value}','fontsize',fs)


