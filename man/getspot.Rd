\name{getspot}

\alias{getanno}
\alias{getspot}
\alias{getsirna}

\title{Mitocheck annotation}

\description{
  Functions to access Mitocheck screen annotation.
}

\usage{
getanno(spot=NULL, sirna=NULL, field="hgnc")
getspot(sirna=NULL, ann=NULL, field="hgnc")
getsirna(spot=NULL, ann=NULL, field="hgnc")
}

\arguments{
  \item{spot}{A numeric indicating a spot ID, synchronised with the global
      Mitocheck \code{tab} object.}
    
  \item{sirna}{A character string indicating a Mitocheck siRNA ID, as 
      referenced in http://www.mitocheck.org.}

  \item{ann}{A character string indicating an annoated value of a
      siRNA ID, from the \code{field} value of the global Mitocheck \code{anno} object.}

  \item{field}{A character vector indicating siRNA annotation field values.
      Possible values include: \code{ensembl}, \code{entrez}, \code{hgnc} and
      \code{getname}. Default is \code{hgnc}.}
}

\value{
  A numeric or character vector.
}

\seealso{
  \code{\link{readspot}}, \code{\link{plotspot}}
}

\author{
  Gregoire Pau, \email{pau.gregoire@gene.com}, 2012
}

\examples{
  ## which siRNAs are targeting the gene CDH1?
  getsirna(ann="CDH1")

  ## which spots contains the siRNA MCO_0026105?
  getspot(sirna="MCO_0026105")

  ## which spots target the gene VIM?
  getspot(ann="VIM")

  ## which gene is targeted by spot 12345?
  getanno(spot=12345)
  getanno(spot=12345, field=c("hgnc", "entrez", "genename"))
}
