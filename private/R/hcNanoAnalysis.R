library(ggplot2)

d <- read.table("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/hcNano/runtimes.txt", header=T)
d$iteration <- as.factor(d$iteration)

p <- ggplot(subset(d, iteration==3), aes(x=nct, y=runtime.min))
p <- p + geom_point() + geom_line() 
p <- p + geom_smooth(se=F, linetype="dashed")
p <- p + ggtitle("HaplotypeCaller performance as a function of NCT")
print(p)

p <- p + coord_trans(x="log2", y="log2")
print(p)