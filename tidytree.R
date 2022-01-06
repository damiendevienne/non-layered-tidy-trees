# Implementation of the non-layered tidy tree algo of Atze van der Ploeg (2014)
# SOFTWARE – PRACTICE AND EXPERIENCE
# Softw. Pract. Exper. 2014; 44:1467–1484
# Published online 19 July 2013 in Wiley Online Library (wileyonlinelibrary.com). DOI: 10.1002/spe.2213


# x is a tree (phylo)
# 
require(ape)


x<-read.tree(text="(((t1:0.02238364564,t2:0.575968103):0.3,((t7:0.2118051811,t8:0.2788739109):0.0458753889,t5:0.8623356319):0.7917935827,(t6:0.2282042881,(t3:0.5514625001,t2:0.7040172669):0.2897282881):0.3026764626):0.7412839716,((t9:0.3830648817,(t1:0.8488779333,t10:0.8147486809):0.542532298):0.5484476415,t4:0.7387184096):0.7571798237);")

xy<-as.data.frame(plotPhyloCoor(x))
edge<-x$edge


nbt<-Ntip(x)
nbn<-Nnode(x)
maxn<-nbt+nbn

oedge<-edge[match(1:maxn,edge[,2]),1]
segofnodes<-data.frame(x1=xy$xx[oedge], y1=xy$yy, x2=xy$xx,y2=xy$yy)

nodes<-maxn:(nbt+1)
N<-list()
for (n in nodes) {
	N[[n]]<-list()
	childs<-edge[edge[,1]==n,2]
	# childbottom<-childs[order(xy$yy[childs])][1]
	# childtop<-childs[order(xy$yy[childs], decreasing=T)][1]
	desc<-c(childs,unlist(lapply(N[childs], function(x) x$desc)))
	# descbottom<-c(childbottom,unlist(lapply(N[childbottom], function(x) x$descbottom)))
	# desctop<-c(childtop,unlist(lapply(N[childtop], function(x) x$desctop)))

	descseg<-segofnodes[childs,]
	segtop.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segtop)))
	segtop<-getContourFromSegments(segtop.pre, "top")

	segbottom.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segbottom)))
	segbottom<-getContourFromSegments(segbottom.pre, "bottom")
	N[[n]]$desc<-desc
	N[[n]]$segtop<-segtop
	N[[n]]$segbottom<-segbottom	
	# N[[n]]$descbottom<-descbottom
	# N[[n]]$desctop<-desctop
}







getContourFromSegments<-function(seg, which) {
	allx<-sort(unique(c(seg$x1, seg$x2)))
	newx2<-allx[2:length(allx)]
	if (which=="top")	newy2<-	sapply(newx2, function(cx,se) max(se[cx>se$x1&cx<=se$x2,]$y1), se=seg) # cx = current x; se = seg
	if (which=="bottom") newy2<-	sapply(newx2, function(cx,se) min(se[cx>se$x1&cx<=se$x2,]$y1), se=seg) # cx = current x; se = seg
	newx1<-allx[1:(length(allx)-1)]
	newy1<-newy2
	newseg<-data.frame(x1=newx1, y1=newy1, x2=newx2, y2=newy2)
	##TODO: simplify new seg to remove segments "bout à bout"
	return(newseg)
}
