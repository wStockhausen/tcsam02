library(ggplot2)
jitterIt<-function(i=50, l=0, u=100, jitFrac=0.6, n=1000,buffer=0.001){
    lo = l;
    uo = u;
    d = u - l;
    l = l+buffer*d;#shrink lower bound
    u = u-buffer*d;#shrink upper bound
    d = u - l;    #shrink interval
    lp = i - 0.5*d*jitFrac;
    up = i + 0.5*d*jitFrac;
    rp = i + (runif(n)-0.5)*d*jitFrac;
    idx<-rp>u; rp[idx] = lp - (rp[idx]-u);
    idx<-rp<l; rp[idx] = up + (l-rp[idx]);

    dfr<-data.frame(r=rp);

    p<-ggplot(data=dfr,mapping=aes_string(x="r")) + geom_histogram();
    p <- p + coord_cartesian(xlim=c(min(lo,lp),max(uo,up)));
    p <- p + geom_vline(xintercept=lo,colour="red") + geom_vline(xintercept=uo,colour="red");
    p <- p + geom_vline(xintercept=lp,linetype=2) + geom_vline(xintercept=up,linetype=2);
    p <- p + geom_vline(xintercept=i,color="blue")
    print(p)
    return(rp);
}

r<-jitterIt(i=80);

