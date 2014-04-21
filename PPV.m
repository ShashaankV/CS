%(E)xpected
%(O)bserved
%nz=nonzero
%z=zero
function PPV=PPV(bool_nzE,bool_nzO)
    P=sum(bool_nzO); %number of called positives
    TP=sum(bool_nzE.*bool_nzO);
    PPV=TP/P;
end