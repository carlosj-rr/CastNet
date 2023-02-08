par(mfrow=c(1,1),mar=c(4,4,4,4))

files_in<-list.files(pattern = "fitness")

tables_list<-list()

for (i in 1:length(files_in)) {
    tables_list[[i]]<-read.table(files_in[i],sep=",",header=F)
}

ancs<-list(tables_list[[3]],tables_list[[6]])

tips<-list(tables_list[[1]],tables_list[[2]],tables_list[[4]],tables_list[[5]])

plot(1:101,ancs[[1]]$V1,type="l",col="red",ylab="Population mean fitness",xlab="Generation num (X100)",xlim=c(0,200),ylim=c(0.5,1),lwd=3)

lines(1:101,ancs[[2]]$V1, col="green", lwd=2) # lines that are plotted above other lines should be thinner, to show the other one behind

lines(102:202,tips[[1]]$V1, col="purple",lwd=3)

lines(102:202,tips[[2]]$V1, col="cyan",lwd=2)

lines(102:202,tips[[3]]$V1, col="magenta",lwd=3)

lines(102:202,tips[[4]]$V1, col="black",lwd=2)

abline(v=101,col="orange",lwd=1) # a vertical line showing the split in lineages

legend(125,0.790,legend=c("(A,B) ancestor","(C,D) ancestor", "Lineage A","Lineage B","Lineage C","Lineage D","Split"),fill=c("red","green","purple","cyan","magenta","black","orange"),cex=0.7)