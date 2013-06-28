loadMitocheck <- function() {
  if (!exists("tab")) {
    ## load assay layout
    tab <- get(load(system.file("data/tab.rda", package="mitoODEdata")))
    assign("tab", tab, 1)
    
    ## load assay annotation
    tab <- get(load(system.file("data/anno.rda", package="mitoODEdata")))
    assign("anno", anno, 1)
    
    ## load assay start/stop time
    assign("g.tstart", 18, 1)
    assign("g.tstop", g.tstart+48, 1)
  }
}

loadPheno <- function() {
  if (!exists("pheno")) {
    tab <- get(load(system.file("data/pheno.rda", package="mitoODEdata")))
    assign("pheno", pheno, 1)
  }
}

readspot <- function(spot) {
  loadMitocheck()
  prs <- c(plate=tab$plate[spot], replicate=tab$replicate[spot],
           spot=tab$spot[spot])
  plate <- prs["plate"]
  replicate <- prs["replicate"]
  spot <- prs["spot"]
  y <- NULL

  ## cache, to speed up readspot
  if (exists("cache.plate")) {
    if (plate==cache.plate.pr[1] && replicate==cache.plate.pr[2]) y <- cache.plate[[spot]]
  }

  ## if not cached, read spot
  if (is.null(y)) {
    platename <- sprintf("LT%04d_%02d.rda", plate, replicate)
    platepath <- file.path(system.file("data", package="mitoODEdata"), platename)
    cdat <- get(load(platepath))
    cache.plate <- cdat
    cache.plate.pr <- prs[1:2]
    assign('cache.plate',cache.plate,1)
    assign('cache.plate.pr',cache.plate.pr,1)
    y <- cache.plate[[spot]]
  }

  y
}

plotspot <- function(spot) {
  loadMitocheck()
  
  ## constants
  lwd <- 1
  xlab <- "Time after cell seeding (h)"
  ylab <- "Number of cells"
  legend <- "topleft"
  cex <- 1

  ## plot spot
  y <- try(readspot(spot), silent=TRUE)
  if (class(y)=='try-error') cat('warning: cannot read spot=', spot, '\n')
  else {
    nt <- nrow(y)
    vt <- seq(g.tstart, len=nt, by=0.5)
    
    cola <- c("#00000077", "#FF000077", "#00CD0077", "#0000FF77")
    colb <- c("#000000", "#FF0000", "#00CD00", "#0000FF")
    matplot(vt, y, type="l", xlab=xlab, ylab=ylab, lty=1, lwd=lwd, cex=cex, col=colb)
    if (!is.null(legend)) legend(legend,
                                 legend=c('interphase', 'mitotic', 'polynucleated', 'dead'),
                                 col=c(1,2,3,4), lty=1, lwd=lwd, cex=cex, seg.len=1,
                                 bg="#ffffff")
  }
}

getanno <- function(spot=NULL, sirna=NULL, field="hgnc") {
  loadMitocheck()
  if (!is.null(spot)) getanno(sirna=getsirna(spot), field=field)
  else anno[match(sirna, anno$sirna), field]
}

getspot <- function(sirna=NULL, ann=NULL, field="hgnc") {
  loadMitocheck()
  if (!is.null(ann)) getspot(sirna=getsirna(ann=ann))
  else {
    z <- which(tab$sirna %in% sirna)
    unlist(split(z, tab$sirna[z])[sirna], use.names=FALSE)
  }
}

getsirna <- function(spot=NULL, ann=NULL, field="hgnc") {
  loadMitocheck()
  if (!is.null(spot)) tab$sirna[spot]
  else {
    z <- which(anno[,field]%in%ann)
    unlist(split(anno$sirna[z], anno[z, field])[ann], use.names=FALSE)
  }
}
