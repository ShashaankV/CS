function P=lse(y,X)
    b=(inv(X'*X))*(X'*y);
    ss_yy=sum((y-mean(y)).^2);
    mu_X=mean(X);
    [m,n]=size(X);
    ss_xx=sum((X-repmat(mu_X,[m,1])).^2);
    ss_xy=sum((X-repmat(mu_X,[m,1])).*repmat((y-mean(y)),[1,n]));
    s=sqrt((ss_yy-b'.*ss_xy)./(m-2));
    se=s./sqrt(ss_xx);
    df=m-1;
    t=b/se;
    P=2*(tcdf(-abs(t),df));
end