# dada2:::learnErrors

learnErrors(filts, nreads = nreads.learn, multithread = multithread, 
    HOMOPOLYMER_GAP_PENALTY = HOMOPOLYMER_GAP_PENALTY, BAND_SIZE = BAND_SIZE)
{
    if (!is.null(nreads)) {
        warning("The nreads parameter is DEPRECATED. Please update your code with the nbases parameter.")
    }
    NBASES <- 0
    NREADS <- 0
    if (is(fls, "derep")) {
        fls <- list(fls)
    }
    if (is.character(fls) && length(fls) == 1 && dir.exists(fls)) {
        fls <- parseFastqDirectory(fls)
    }
    drps <- vector("list", length(fls))
    if (randomize) {
        fls <- sample(fls)
    }
    for (i in seq_along(fls)) {
        if (is.list.of(fls, "derep")) {
            drps[[i]] <- fls[[i]]
        }
        else {
            drps[[i]] <- derepFastq(fls[[i]], qualityType = qualityType)
        }
        NREADS <- NREADS + sum(drps[[i]]$uniques)
        NBASES <- NBASES + sum(drps[[i]]$uniques * nchar(names(drps[[i]]$uniques)))
        if (is.null(nreads) && NBASES > nbases) {
            break
        }
        if (!is.null(nreads) && NREADS > nreads) {
            break
        }
    }
    drps <- drps[1:i]
    if (is.logical(verbose) || verbose > 0) {
        cat(NBASES, "total bases in", NREADS, "reads from", i, 
            "samples will be used for learning the error rates.\n")
    }
    dds <- dada(drps, err = NULL, errorEstimationFunction = errorEstimationFunction, 
        selfConsist = TRUE, multithread = multithread, verbose = verbose, 
        MAX_CONSIST = MAX_CONSIST, OMEGA_C = OMEGA_C, ...)
    return(getErrors(dds, detailed = TRUE))
}