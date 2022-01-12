# Implementation of the non-layered tidy tree algo of Atze van der Ploeg (2014)
# SOFTWARE – PRACTICE AND EXPERIENCE
# Softw. Pract. Exper. 2014; 44:1467–1484
# Published online 19 July 2013 in Wiley Online Library (wileyonlinelibrary.com). DOI: 10.1002/spe.2213

# This version is a simplistic one, that does not 
# allow for different orientations, colors, lty, lwd, etc. 
# for full options, use type="tidy" in plot.phylo() 
# in the "ape" package (Paradis et al.)


# x is a tree (phylo)
# 
require(ape)



plot.tidy<-function(x, show.tip.label=T, cex=1, x.lim = NULL, y.lim = NULL, verbose=T, ...) {
	xy<-as.data.frame(plotPhyloCoor(x))
	#get basic tree info
	xe<-x$edge
	Ntip<-Ntip(x)
	Nnode<-Nnode(x)
	xx<-xy$xx
	yy<-xy$yy

    getStringLengthbyTip <- function(x,lab,sin,cex) {
        s <- strwidth(lab, "inches", cex = cex)
        lim<-getLimit(x,lab,sin,cex)
        alp<-lim/sin
        snew<-s*alp
        return(snew)
    }
    ### END TIDY

    ## Function to compute the axis limit
    ## x: vector of coordinates, must be positive (or at least the largest value)
    ## lab: vector of labels, length(x) == length(lab)
    ## sin: size of the device in inches
    getLimit <- function(x, lab, sin, cex) {
        s <- strwidth(lab, "inches", cex = cex) # width of the tip labels
        ## if at least one string is larger than the device,
        ## give 1/3 of the plot for the tip labels:
        if (any(s > sin)) return(1.5 * max(x))
        Limit <- 0
        while (any(x > Limit)) {
            i <- which.max(x)
            ## 'alp' is the conversion coeff from inches to user coordinates:
            alp <- x[i]/(sin - s[i])
            Limit <- x[i] + alp*s[i]
            x <- x + alp*s
        }
        Limit
    }

	if (!show.tip.label) {
		yy<-tidy.xy(xe, Ntip, Nnode, xx, yy, verbose)
	} else { #we add to xx the size taken by labels, so that tidying considers labels 
		xx.tips<-xx[1:Ntip]
	    pin1 <- par("pin")[1] # width of the device in inches               
	    lab.strlength<-getStringLengthbyTip(xx.tips,x$tip.label, pin1, cex) #size of lab strings
	    xx2<-xx
	    xx2[1:Ntip]<-xx2[1:Ntip]+lab.strlength
	    yy<-tidy.xy(xe, Ntip, Nnode, xx2, yy, verbose) #compress taking labels into account
	}
	xy<-data.frame(xx=xx, yy=yy)
	if (is.null(x.lim)) {
        xx.tips <- xx[1:Ntip]# * 1.04
        if (show.tip.label) {
            pin1 <- par("pin")[1] # width of the device in inches
            tmp <- getLimit(xx.tips, x$tip.label, pin1, cex)
        } else tmp <- max(xx.tips)
        x.lim <- c(0, tmp)
	}
	if (is.null(y.lim)) {
		y.lim <- c(1, max(yy[1:Ntip]))
	}
	plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "",
                 ylab = "", axes = FALSE, asp = NA, ...)

	plottree(xe,xy)
	if (show.tip.label) {
		text(xx[1:Ntip], yy[1:Ntip], x$tip.label, cex=cex, pos=4, offset=0)
	}

}

segplot<-function(mat, col="red", lwd=2) {
	segments(mat[,1],mat[,2], mat[,3], mat[,4], col=col, lwd=lwd)
}


plottree<-function(edge, xy) {
#	plot(xy, type="n", frame.plot=FALSE, axes=F, xlab="", ylab="", xlim=xlim, ylim=ylim)	
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


tidy.xy<-function(edge, Ntip, Nnode, xx, yy, verbose=F) { 

    yynew<-yy #will be updated to get the new y coordinates after tidying

    initialrange<-diff(range(yy)) #for computing compression. Remove ? 

    oedge<-edge[match(seq_len(Ntip+Nnode),edge[,2]),1] # ordered edges 
    segofnodes<-data.frame(x1=xx[oedge], y1=yy, x2=xx,y2=yy) # segment associated to each node

    nodes<-order(xx, decreasing=T)

    GetContourPairsFromSegments<-function(seg, which) {
        allx<-sort(unique(c(seg$x1, seg$x2)))
        newx2<-allx[2:length(allx)]
        if (which=="top")   {
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

    N<-list() ##will contain all info for each node.
    for (n in nodes) {
        N[[n]]<-list()
        childs<-edge[edge[,1]==n,2]
        childs.ord<-childs[order(yy[childs])] #childs orderd by y values
        desc<-c(childs,unlist(lapply(N[childs], function(x) x$desc)))

        diffiny<-0

        if (n>Ntip) { #we are in a node
            oldyofcurrentnode<-yynew[n]
            for (nn in 2:length(childs.ord)) {
                top<-N[[childs.ord[nn-1]]]$segtop
                bot<-N[[childs.ord[nn]]]$segbottom
                mindist<-GetMinDistBetweenContours(top, bot)
                if (mindist != 1) { ## There is room for tidying or untidy up if branches are tangled 
                    mod<-mindist-1
                    N[[childs.ord[nn]]]$segbottom[,c(2,4)]<-N[[childs.ord[nn]]]$segbottom[,c(2,4)]-mod
                    N[[childs.ord[nn]]]$segtop[,c(2,4)]<-N[[childs.ord[nn]]]$segtop[,c(2,4)]-mod
                
                    yynew[c(childs.ord[nn], N[[childs.ord[nn]]]$desc)]<-yynew[c(childs.ord[nn], N[[childs.ord[nn]]]$desc)]-mod
                }
            }
            newyofcurrentnode<-mean(range(yynew[childs.ord]))
            yynew[n]<-newyofcurrentnode
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
    yynew<-yynew-(min(yynew)-1) ## so that min(y)=1 always
    finalrange<-diff(range(yynew)) #for computing compression. Remove ? 
    compression<-((initialrange-finalrange)/initialrange)*100 # Remove ? 
    if (verbose) print(paste("Compression: ", round(compression,2),"%", sep="")) #Remove ?
    return(yynew)
    #       
}
