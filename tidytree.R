# Implementation of the non-layered tidy tree algo of Atze van der Ploeg (2014)
# SOFTWARE – PRACTICE AND EXPERIENCE
# Softw. Pract. Exper. 2014; 44:1467–1484
# Published online 19 July 2013 in Wiley Online Library (wileyonlinelibrary.com). DOI: 10.1002/spe.2213


# x is a tree (phylo)
# 
require(ape)


x<-read.tree(text="(((t1:0.02238364564,t2:0.575968103):0.3,((t7:0.2118051811,t8:0.2788739109):0.0458753889,t5:0.8623356319):0.7917935827,(t6:0.2282042881,(t3:0.5514625001,t2:0.7040172669):0.2897282881):0.3026764626):0.7412839716,((t9:0.3830648817,(t1:0.8488779333,t10:0.8147486809):0.542532298):0.5484476415,t4:0.7387184096):0.7571798237);")



plot.tidy(x) {

	xy<-as.data.frame(plotPhyloCoor(x))
	
	########## FOR DEV PURPOSE ONLY
	plot(x)
	nodelabels()
	tiplabels()
	####################

	#get basic tree info
	edge<-x$edge
	nbt<-Ntip(x)
	nbn<-Nnode(x)
	maxn<-nbt+nbn

	#prepare data
	oedge<-edge[match(1:maxn,edge[,2]),1]
	segofnodes<-data.frame(x1=xy$xx[oedge], y1=xy$yy, x2=xy$xx,y2=xy$yy)
	
	nodes<-order(xy$xx, decreasing=T)
	#first traversal: get contours for each node
	N<-list()
	for (n in nodes) {
		N[[n]]<-list()
		childs<-edge[edge[,1]==n,2]
		childs.ord<-childs[order(xy$yy[childs])] #childs orderd by y values (for 2nd traversal)
		# childbottom<-childs[order(xy$yy[childs])][1]
		# childtop<-childs[order(xy$yy[childs], decreasing=T)][1]
		desc<-c(childs,unlist(lapply(N[childs], function(x) x$desc)))
		# descbottom<-c(childbottom,unlist(lapply(N[childbottom], function(x) x$descbottom)))
		# desctop<-c(childtop,unlist(lapply(N[childtop], function(x) x$desctop)))

		descseg<-segofnodes[childs,]
		##we add to descseg the seg leading to the current node 
		descseg<-rbind(descseg, segofnodes[n,])

		segtop.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segtop)))
		segtop<-getContourFromSegments(segtop.pre, "top")

		segbottom.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segbottom)))
		segbottom<-getContourFromSegments(segbottom.pre, "bottom")

		N[[n]]$childs<-childs.ord
		N[[n]]$desc<-desc
		N[[n]]$segtop<-segtop
		N[[n]]$segbottom<-segbottom	



		## at each step we can look at the childs (they have already been treated by def)
		## et compute the movement (mod)
		if (n>nbt) { #we are in a node
			#for each child c in the correct order (childs.ord)
			#compare the top contour of c the bottom contour of c+1
			for (nn in 2:length(childs.ord)) {
				top<-N[[childs.ord[nn-1]]]$segtop
				bot<-N[[childs.ord[nn]]]$segbottom
				segplot(top, col="red")
				segplot(bot, col="green")
#				compare(N[[childs.ord[nn-1]]]$segtop, N[[childs.ord[nn]]]$segbottom)
				scan()
				plot(x)
				nodelabels()
				tiplabels()

			}
		}
	}
	# second traversal: compute for each subtree it's possible movement, given 
	# the distance between its bottom contour and the top contour of its bottom sibling subtree
	# at each node how much it should move down. 
	# We explore the tree by nodes in decreasing x values. We always move nodes down.
	# Childs are taken in in creasing y value order. (Because we always move nodes down)
}


d<-NULL
topi<-1
boti<-1
while((topi<=nrow(top))&(boti<=nrow(bot))) {
	d<-c(d,bot[boti,]$y1-top[topi,]$y1)
	if (bot[boti,]$x2<top[topi,]$x2) boti<-boti+1
	else topi<-topi+1
	print(topi)
	print(boti)
}



getContourFromSegments<-function(seg, which) {
	allx<-sort(unique(c(seg$x1, seg$x2)))
	newx2<-allx[2:length(allx)]
	if (which=="top")	{
		newy2i<-sapply(newx2, function(cx,se) which(cx>se$x1&cx<=se$x2)[which.max(se$y1[which(cx>se$x1&cx<=se$x2)])], se=seg)
	}
	if (which=="bottom") {
		newy2i<-sapply(newx2, function(cx,se) which(cx>se$x1&cx<=se$x2)[which.min(se$y1[which(cx>se$x1&cx<=se$x2)])], se=seg)
	}

	newx1<-allx[1:(length(allx)-1)]
	newy1i<-newy2i
	#we simplify segments by merging thsoe on same horiz  (bout à bout)
	where2mergei<-which((newy1i[2:length(newy1i)]-newy2i[1:(length(newy1i)-1)])==0)
	if(length(where2mergei)>0) {
		newx1<-newx1[-(where2mergei+1)]
		newy1i<-newy1i[-(where2mergei+1)]
		newx2<-newx2[-(where2mergei)]
		newy2i<-newy2i[-(where2mergei)]
	}

	newy1ok<-seg$y1[newy1i]
	newy2ok<-newy1ok
	newseg<-data.frame(x1=newx1, y1=newy1ok, x2=newx2, y2=newy2ok)
	##TODO: simplify new seg to remove segments "bout à bout"
	return(newseg)
}


segplot<-function(mat, col="red", lwd=2) {
	segments(mat[,1],mat[,2], mat[,3], mat[,4], col=col, lwd=lwd)
}