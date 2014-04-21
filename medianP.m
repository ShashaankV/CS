function mu=medianP(xhat,G,miss_value,y)
    A=zscore_sv(G,miss_value,'zero');
    indxhat=find(xhat~=0);
    P=[];
    for i=1:length(indxhat)
        P(i)=lse(y,A(:,indxhat(i)));
    end
    mu=median(P);
end