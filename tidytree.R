# Implementation of the non-layered tidy tree algo of Atze van der Ploeg (2014)
# SOFTWARE – PRACTICE AND EXPERIENCE
# Softw. Pract. Exper. 2014; 44:1467–1484
# Published online 19 July 2013 in Wiley Online Library (wileyonlinelibrary.com). DOI: 10.1002/spe.2213


# x is a tree (phylo)
# 
require(ape)


GetMinDistBetweenContours<- function(topcontour, bottomcontour) {
	## efficient way to compare top and bottom contour by only looking
	## at necessary pairs (see original publiction)
	d<-NULL
	topi<-1
	boti<-1
	while((topi<=nrow(topcontour))&(boti<=nrow(bottomcontour))) {
		d<-c(d,bottomcontour[boti,]$y1-topcontour[topi,]$y1)
		if (bottomcontour[boti,]$x2<topcontour[topi,]$x2) boti<-boti+1
		else topi<-topi+1
	}
	return(min(d))
}


GetContourPairsFromSegments<-function(seg, which) {
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


plottree<-function(edge, xy, xlim, ylim) {
	plot(xy, type="n", frame.plot=FALSE, axes=F, xlab="", ylab="", xlim=xlim, ylim=ylim)	
	X1<-xy$xx[edge[,1]]
	Y1<-xy$yy[edge[,2]]
	X2<-xy$xx[edge[,2]]
	Y2<-xy$yy[edge[,2]]
	segments(X1,Y1,X2,Y2)
	XX1<-X1
	YY1<-xy$yy[edge[,1]]
	XX2<-XX1
	YY2<-Y2
	segments(XX1, YY1,XX2,YY2)
}



silent<-function(x) {

#	par(mfrow=c(1,2))


#	x<-read.tree(text="(((t1:0.02238364564,t2:0.575968103):0.3,((t7:0.2118051811,t8:0.2788739109):0.0458753889,t5:0.8623356319):0.7917935827,(t6:0.2282042881,(t3:0.5514625001,t2:0.7040172669):0.2897282881):0.3026764626):0.7412839716,((t9:0.3830648817,(t1:0.8488779333,t10:0.8147486809):0.542532298):0.5484476415,t4:0.7387184096):0.7571798237);")

	xy<-as.data.frame(plotPhyloCoor(x))
	
	#get basic tree info
	edge<-x$edge


	##PLOT FUN
	x.lim<-range(xy$xx)
	y.lim<-range(xy$yy)

	par(mfrow=c(1,2))
	plottree(edge, xy, xlim=x.lim, ylim=y.lim)



	Ntip<-Ntip(x)
	Nnode<-Nnode(x)
	maxn<-Ntip+Nnode

	#prepare data
	oedge<-edge[match(seq_len(maxn),edge[,2]),1]
	segofnodes<-data.frame(x1=xy$xx[oedge], y1=xy$yy, x2=xy$xx,y2=xy$yy)
	
	nodes<-order(xy$xx, decreasing=T)
	N<-list()
	for (n in nodes) {
		N[[n]]<-list()
		childs<-edge[edge[,1]==n,2]
		childs.ord<-childs[order(xy$yy[childs])] #childs orderd by y values (for 2nd traversal)
		desc<-c(childs,unlist(lapply(N[childs], function(x) x$desc)))

		diffiny<-0

		if (n>Ntip) { #we are in a node
			oldyofcurrentnode<-xynew$yy[n]
			for (nn in 2:length(childs.ord)) {
				miniychild<-NULL
				maxiychild<-NULL
				top<-N[[childs.ord[nn-1]]]$segtop
				bot<-N[[childs.ord[nn]]]$segbottom
				mindist<-GetMinDistBetweenContours(top, bot)
				if (mindist > 1) { ## There is room for tidying
					mod<-mindist-1
					N[[childs.ord[nn]]]$segbottom[,c(2,4)]<-N[[childs.ord[nn]]]$segbottom[,c(2,4)]-mod
					N[[childs.ord[nn]]]$segtop[,c(2,4)]<-N[[childs.ord[nn]]]$segtop[,c(2,4)]-mod
				
					xynew$yy[c(childs.ord[nn], N[[childs.ord[nn]]]$desc)]<-xynew$yy[c(childs.ord[nn], N[[childs.ord[nn]]]$desc)]-mod
				}
			}
			newyofcurrentnode<-mean(range(xynew$yy[childs.ord]))
			xynew$yy[n]<-newyofcurrentnode
			diffiny<-oldyofcurrentnode-newyofcurrentnode
		}

		descseg<-segofnodes[n,]
		descseg[,c(2,4)]<-descseg[,c(2,4)]-diffiny
		segtop.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segtop)))
		segtop<-GetContourPairsFromSegments(segtop.pre, "top")
		segbottom.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segbottom)))
		segbottom<-GetContourPairsFromSegments(segbottom.pre, "bottom")

		N[[n]]$childs<-childs.ord
		N[[n]]$desc<-desc
		N[[n]]$segtop<-segtop
		N[[n]]$segbottom<-segbottom	

	} 
plottree(edge, xynew, xlim=x.lim, ylim=y.lim)
}



tidy.xy<-function(edge, Ntip, Nnode, xx, yy) {
	xynew<-data.frame(xx=xx,yy=yy) #will be updated to get the new coordinates at the end

	oedge<-edge[match(seq_len(Ntip+Nnode),edge[,2]),1] # ordered edges 
	segofnodes<-data.frame(x1=xx[oedge], y1=yy, x2=xx,y2=yy) # segment associated to each node

	nodes<-order(xx, decreasing=T)
	N<-list() ##will contain all info for each node.

	for (n in nodes) {
		N[[n]]<-list()
		childs<-edge[edge[,1]==n,2]
		childs.ord<-childs[order(yy[childs])] #childs orderd by y values (for 2nd traversal)
		desc<-c(childs,unlist(lapply(N[childs], function(x) x$desc)))

		diffiny<-0

		if (n>Ntip) { #we are in a node
			oldyofcurrentnode<-xynew$yy[n]
			for (nn in 2:length(childs.ord)) {
				miniychild<-NULL
				maxiychild<-NULL
				top<-N[[childs.ord[nn-1]]]$segtop
				bot<-N[[childs.ord[nn]]]$segbottom
				mindist<-GetMinDistBetweenContours(top, bot)
				if (mindist > 1) { ## There is room for tidying
					mod<-mindist-1
					N[[childs.ord[nn]]]$segbottom[,c(2,4)]<-N[[childs.ord[nn]]]$segbottom[,c(2,4)]-mod
					N[[childs.ord[nn]]]$segtop[,c(2,4)]<-N[[childs.ord[nn]]]$segtop[,c(2,4)]-mod
				
					xynew$yy[c(childs.ord[nn], N[[childs.ord[nn]]]$desc)]<-xynew$yy[c(childs.ord[nn], N[[childs.ord[nn]]]$desc)]-mod
				}
			}
			newyofcurrentnode<-mean(range(xynew$yy[childs.ord]))
			xynew$yy[n]<-newyofcurrentnode
			diffiny<-oldyofcurrentnode-newyofcurrentnode
		}

		descseg<-segofnodes[n,]
		descseg[,c(2,4)]<-descseg[,c(2,4)]-diffiny
		segtop.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segtop)))
		segtop<-GetContourPairsFromSegments(segtop.pre, "top")
		segbottom.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segbottom)))
		segbottom<-GetContourPairsFromSegments(segbottom.pre, "bottom")

		N[[n]]$childs<-childs.ord
		N[[n]]$desc<-desc
		N[[n]]$segtop<-segtop
		N[[n]]$segbottom<-segbottom	

	}
	return(xynew)		
}




tidy.plot <- function(edge, Ntip, Nnode, xx, yy, horizontal,
                           edge.color = NULL, edge.width = NULL,
                           edge.lty = NULL,
                           node.color = NULL, node.width = NULL,
                           node.lty = NULL) 
{

	xynew<-data.frame(xx=xx,yy=yy) #will be updated to get the new coordinates at the end

	oedge<-edge[match(seq_len(Ntip+Nnode),edge[,2]),1] #ordered edges 
	segofnodes<-data.frame(x1=xx[oedge], y1=yy, x2=xx,y2=yy) 

	nodes<-order(xx, decreasing=T)
	N<-list() ##will contain all info for each node.

	for (n in nodes) {
		N[[n]]<-list()
		childs<-edge[edge[,1]==n,2]
		childs.ord<-childs[order(yy[childs])] #childs orderd by y values (for 2nd traversal)
		desc<-c(childs,unlist(lapply(N[childs], function(x) x$desc)))

		diffiny<-0

		if (n>Ntip) { #we are in a node
			oldyofcurrentnode<-xynew$yy[n]
			for (nn in 2:length(childs.ord)) {
				miniychild<-NULL
				maxiychild<-NULL
				top<-N[[childs.ord[nn-1]]]$segtop
				bot<-N[[childs.ord[nn]]]$segbottom
				mindist<-GetMinDistBetweenContours(top, bot)
				if (mindist > 1) { ## There is room for tidying
					mod<-mindist-1
					N[[childs.ord[nn]]]$segbottom[,c(2,4)]<-N[[childs.ord[nn]]]$segbottom[,c(2,4)]-mod
					N[[childs.ord[nn]]]$segtop[,c(2,4)]<-N[[childs.ord[nn]]]$segtop[,c(2,4)]-mod
				
					xynew$yy[c(childs.ord[nn], N[[childs.ord[nn]]]$desc)]<-xynew$yy[c(childs.ord[nn], N[[childs.ord[nn]]]$desc)]-mod
				}
			}
			newyofcurrentnode<-mean(range(xynew$yy[childs.ord]))
			xynew$yy[n]<-newyofcurrentnode
			diffiny<-oldyofcurrentnode-newyofcurrentnode
		}

		descseg<-segofnodes[n,]
		descseg[,c(2,4)]<-descseg[,c(2,4)]-diffiny
		segtop.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segtop)))
		segtop<-GetContourPairsFromSegments(segtop.pre, "top")
		segbottom.pre<-rbind(descseg, do.call(rbind, lapply(N[childs], function(x) x$segbottom)))
		segbottom<-GetContourPairsFromSegments(segbottom.pre, "bottom")

		N[[n]]$childs<-childs.ord
		N[[n]]$desc<-desc
		N[[n]]$segtop<-segtop
		N[[n]]$segbottom<-segbottom	

	} 

#	plot(xy, type="n", frame.plot=FALSE, axes=F, xlab="", ylab="", xlim=xlim, ylim=ylim)	
	#horizontal segments
	X1<-xynew$xx[edge[,1]]
	Y1<-xynew$yy[edge[,2]]
	X2<-xynew$xx[edge[,2]]
	Y2<-xynew$yy[edge[,2]]
	segments(X1,Y1,X2,Y2)
	##vertical segments
	XX1<-X1
	YY1<-xynew$yy[edge[,1]]
	XX2<-XX1
	YY2<-Y2
	segments(XX1, YY1,XX2,YY2)
}







# X<-seq(100,3000, by=100)
# trees<-sapply(X, function(x) rtree(x), simplify=F)
# Y<-unlist(lapply(trees, function(x) system.time(silent(x))[[1]]))
# plot(X,Y, type="o")
