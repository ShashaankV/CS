function maf=calc_maf(G,miss_value)
    [m,n]=size(G);
    for i=1:n
        ind=G(:,i)~=miss_value;
        maf(i)=sum(G(ind,i))/(2*sum(ind)); %ind is logicial, maximum number of minor allelese is twice the number of subjects since biallelic
    end
end