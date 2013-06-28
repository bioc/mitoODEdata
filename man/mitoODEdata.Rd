\name{mitoODEdata}

\alias{mitoODEdata}
\alias{cdat}
\alias{anno}
\alias{pheno}
\alias{tab}
\alias{readspot}	
\alias{plotspot}

\title{Access to the Mitocheck assay data}

\description{
 
The Mitocheck screen \cite{mitocheck} is a time-lapse imaging assay that employed
small-interfering RNAs (siRNAs) to test
the implication of human genes in transient biological processes such as cell division or
migration genome-wide. In this experiment, HeLa cells stably expressing
core histone 2B tagged with green fluorescent protein (GFP) were seeded on siRNA-spotted
slides, incubated for 18 h and imaged with automated fluorescence microscopy for
48~h. Video sequences of cell populations on each siRNA-spot were analysed by image
segmentation, and at each frame, each individual cell was categorised into one of 16
morphological classes mostly related to cell division.

The \code{mitoODE} package implements a modelling by differential
equations of cellular populations \cite{mitoODE}, to quantify the
phenotypic effect induced by siRNA treatments in the Mitocheck
screen. The package includes the code to fit any time course data to
the model and the scripts used to generate the figures and results
presented in the paper.

The \code{mitoODEdata} package, the experimental companion package of
\code{mitoODE}, contains the screen data and methods to
access the Mitocheck assay layout, siRNA annotation, time-lapse cell
counts and the fitted phenotypes for each spot. Four cell types are considered: interphase (referred in the
Mitocheck paper as: Interphase, Large, Elongated, Folded, Hole,
SmallIrregular or UndefinedCondensed), mitotic (Metaphase, Anaphase,
MetaphaseAlignment, Prometaphase or ADCCM), polynucleated (Shape1,
Shape3, Grape) and apoptotic (Apoptosis).

}

\usage{
readspot(spot)
plotspot(spot)
}

\arguments{
  \item{spot}{A numeric indicating a spot ID, synchronised with the global
      Mitocheck \code{tab} object.}
}

\value{
  The function \code{readspot} returns a matrix containing the number of
cells of a given type (interphase ``i'', mitotic ``m'', polynucleated ``s'' and
apoptotic ``a'') per frame. The first frame (e.g. row) was acquired 18 h after
cell seeding and following frames were acquired every 30 minutes. 
}

\details{
  Loading the package \code{mitoODEdata} loads the Mitocheck screen
  annotation variables \code{tab} and \code{anno} in the global environement.
  The object \code{tab} is a data.frame containing spot metadata, including: 
  \code{plate} (plate number), \code{replicate} (replicate
  number), \code{spot} (spot number within the plate), \code{qc}
  (original quality control from the paper), \code{type} (spot type) and
  \code{sirna} (spot siRNA ID). The object \code{anno} is a data.frame
  containing the siRNA to gene mapping, including: \code{sirna} (siRNA
  ID), \code{ensembl} (target Ensembl gene ID), \code{hgnc} (target HGNC
  gene symbol), \code{entrez} (target Entrez gene ID), \code{genename}
  (target HGNC gene name).
}

\seealso{
  \code{\link{getspot}}, \code{\link{getsirna}}, \code{\link{getanno}}
}

\references{
Pau G, Walter T, Neumann B,  Heriche JK, Ellenberg J, and Huber W (2013)
{D}ynamical modelling of phenotypes in a genome-wide RNAi
  live-cell imaging assay. Submitted.\\

Neumann B, Walter T, Heriche JK, Bulkescher J, Erfle H, et~al. (2010)
  {P}henotypic profiling of the human genome by time-lapse microscopy reveals
  cell division genes. Nature Nature 464: 721--727.
}


\author{
  Gregoire Pau, \email{pau.gregoire@gene.com}, 2012
}

\examples{
 ## read spot
 spotid <- getspot(ann="FGFR2")[1]
 y <- readspot(spotid)
 y[1:10,]

 ## plot spot
 plotspot(spotid)
}


