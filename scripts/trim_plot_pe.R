library(ggplot2)
library(reshape2)
library(scales)

data <- read.table( "~/Downloads/trim.tab", header=TRUE )
rownames(data) <- data[,1]
data[,1] <- NULL
data <- t(data)

x <- data.frame( Sample=colnames(data),
                Total_Reads=as.numeric(as.matrix(data["TotalReads",])), 
                                Both_Surviving=as.numeric(as.matrix(data["BothSurviving",])),
                                                LeftMate_Only=as.numeric(as.matrix(data["LeftMateOnly",])),
                                                                RightMate_Only=as.numeric(as.matrix(data["RightMateOnly",])),
                                                                                Dropped=as.numeric(as.matrix(data["Dropped",]))
                )
x1 <- melt(x, id.var="Sample")

upper_limit <- max(x$Total_Reads)
limits <- seq( 0, upper_limit, length.out=10)
colors <- c(Total_Reads="Grey", Both_Surviving="Blue", LeftMate_Only="Green",
            RightMate_Only="Yellow", Dropped="Red")

q <- ggplot(x1, aes(x=Sample, y=value, fill=variable)) + geom_bar( stat =
                                                                  "identity",
                                                                  position="dodge")
+ scale_y_continuous("",limits=c(0,upper_limit), labels=comma, breaks=limits) +
scale_fill_manual(values=colors) + labs( title="Trimmomatic Report\n\n", x =
                                        "Sample Names", y="") +
guides(fill=guide_legend(title=NULL)) + theme_bw() 

q + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=10))


