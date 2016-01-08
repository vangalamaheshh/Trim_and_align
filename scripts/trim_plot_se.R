library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs( trailingOnly = TRUE )
data <- read.table( args[1],header=TRUE, row.names=1)

x <- data.frame( Sample=rownames(data),
                Total_Reads=as.numeric(as.matrix(data[,"TotalReads"])), 
                Surviving=as.numeric(as.matrix(data[,"Surviving"])),
                Dropped=as.numeric(as.matrix(data[,"Dropped"]))
                )
x1 <- melt(x, id.var="Sample")

png( args[2], width = 8, height = 8, unit="in",res=300 )
upper_limit <- max(x$Total_Reads)
limits <- seq( 0, upper_limit, length.out=10)
colors <- c(Total_Reads="Grey", Surviving="Blue", Dropped="Red")

q <- ggplot(x1, aes(x=Sample, y=value, fill=variable)) + geom_bar( stat = "identity", position="dodge") 
q + scale_y_continuous("",limits=c(0,upper_limit), labels=comma, breaks=limits) +
scale_fill_manual(values=colors) + labs( title="Trimmomatic Report\n\n", x = "Sample Names", y="") +
guides(fill=guide_legend(title=NULL)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=10))

dev.off()

