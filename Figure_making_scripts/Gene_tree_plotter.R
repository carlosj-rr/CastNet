tree_files<-list.files(pattern = "treefile")

library(ape)

tr_list<-list()

for (i in 1:length(tree_files)) { tr_list[[i]]<-read.tree(tree_files[i])}

windowsFonts(C="Courier New")

names<-c(); for (i in 1:5) {names<-c(names,paste("Gene",i))}

par(mfrow=c(2,3),mar=c(1.0,1.0,1.0,1.0))

for (i in 1:length(tr_list)) {
    plot.phylo(tr_list[[i]],x.lim = c(0,0.5),main=names[i],edge.width = 3,cex=1.2, family="C")
    if (i == 1) {
        add.scale.bar(length=0.1, x=0.08, y=1.4)
    }
}
