# dada2:::derepFastq

derepFastq(fls[[i]], qualityType = qualityType)
{
    if (!is.character(fls)) {
        stop("File paths must be provided in character format.")
    }
    if (length(fls) == 1 && dir.exists(fls)) {
        fls <- parseFastqDirectory(fls)
    }
    if (!all(file.exists(fls))) {
        stop("Not all provided files exist.")
    }
    rval <- list()
    for (i in seq_along(fls)) {
        fl <- fls[[i]]
        if (verbose) {
            message("Dereplicating sequence entries in Fastq file: ", 
                fl, appendLF = TRUE)
        }
        f <- FastqStreamer(fl, n = n)
        suppressWarnings(fq <- yield(f, qualityType = qualityType))
        out <- qtables2(fq)
        derepCounts <- out$uniques
        derepQuals <- out$cum_quals
        derepMap <- out$map
        while (length(suppressWarnings(fq <- yield(f, qualityType = qualityType)))) {
            newniques = alreadySeen = NULL
            if (verbose) {
                message(".", appendLF = FALSE)
            }
            out <- qtables2(fq)
            if (ncol(out$cum_quals) > ncol(derepQuals)) {
                derepQuals <- cbind(derepQuals, matrix(NA, nrow = nrow(derepQuals), 
                  ncol = (ncol(out$cum_quals) - ncol(derepQuals))))
            }
            else if (ncol(out$cum_quals) < ncol(derepQuals)) {
                out$cum_quals <- cbind(out$cum_quals, matrix(NA, 
                  nrow = nrow(out$cum_quals), ncol = (ncol(derepQuals) - 
                    ncol(out$cum_quals))))
            }
            alreadySeen <- names(out$uniques) %in% names(derepCounts)
            if (any(alreadySeen)) {
                sqnms = names(out$uniques)[alreadySeen]
                derepCounts[sqnms] <- derepCounts[sqnms] + out$uniques[sqnms]
                derepQuals[sqnms, ] <- derepQuals[sqnms, ] + 
                  out$cum_quals[sqnms, ]
            }
            if (!all(alreadySeen)) {
                derepCounts <- c(derepCounts, out$uniques[!alreadySeen])
                derepQuals <- rbind(derepQuals, out$cum_quals[!alreadySeen, 
                  , drop = FALSE])
            }
            new2old <- match(names(out$uniques), names(derepCounts))
            if (any(is.na(new2old))) 
                warning("Failed to properly extend uniques.")
            derepMap <- c(derepMap, new2old[out$map])
        }
        derepQuals <- derepQuals/derepCounts
        ord <- order(derepCounts, decreasing = TRUE)
        derepCounts <- derepCounts[ord]
        derepQuals <- derepQuals[ord, , drop = FALSE]
        derepMap <- match(derepMap, ord)
        if (verbose) {
            message("Encountered ", length(derepCounts), " unique sequences from ", 
                sum(derepCounts), " total sequences read.")
        }
        close(f)
        derepO <- list(uniques = derepCounts, quals = derepQuals, 
            map = derepMap)
        derepO <- as(derepO, "derep")
        rval[[i]] <- derepO
    }
    if (length(rval) == 1) {
        rval <- rval[[1]]
    }
    else {
        if (is.null(names(fls))) {
            names(rval) <- basename(fls)
        }
        else {
            names(rval) <- names(fls)
        }
    }
    return(rval)
}