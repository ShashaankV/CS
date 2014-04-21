%(E)xpected
%(O)bserved
%nz=nonzero
%z=zero
function FPR=FPR(bool_nzE,bool_nzO)
    EN=sum(1-bool_nzE); %number of expected negatives
    FP=sum(bool_nzO.*(1-bool_nzE));
    FPR=FP/EN; 
end