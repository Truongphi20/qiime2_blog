# dada2:::removeBimeraDenovo

removeBimeraDenovo(seqtab, method = chimeraMethod, minFoldParentOverAbundance = minParentFold, 
    allowOneOff = allowOneOff, multithread = multithread)
{
    if (is(unqs, "dada") || is(unqs, "derep") || is(unqs, "data.frame")) {
        unqs <- list(unqs)
    }
    if (!is.list(unqs)) {
        unqs <- list(unqs)
    }
    if ("tableMethod" %in% names(list(...))) {
        stop("DEFUNCT: The tableMethod argument has been replaced by the method argument. Please update your code.")
    }
    outs <- list()
    for (i in seq_along(unqs)) {
        if (is.integer(unqs[[i]]) && length(names(unqs[[i]])) != 
            0 && !any(is.na(names(unqs[[i]])))) {
            bim <- isBimeraDenovo(unqs[[i]], ..., verbose = verbose)
            outs[[i]] <- unqs[[i]][!bim]
        }
        else if (is(unqs[[i]], "dada")) {
            bim <- isBimeraDenovo(unqs[[i]], ..., verbose = verbose)
            outs[[i]] <- unqs[[i]]$denoised[!bim]
        }
        else if (is(unqs[[i]], "derep")) {
            bim <- isBimeraDenovo(unqs[[i]], ..., verbose = verbose)
            outs[[i]] <- unqs[[i]]$uniques[!bim]
        }
        else if (is.data.frame(unqs[[i]]) && all(c("sequence", 
            "abundance") %in% colnames(unqs[[i]]))) {
            bim <- isBimeraDenovo(unqs[[i]], ..., verbose = verbose)
            outs[[i]] <- unqs[[i]][!bim, ]
        }
        else if (is.matrix(unqs[[i]]) && !any(is.na(colnames(unqs[[i]])))) {
            if (method == "pooled") {
                bim <- isBimeraDenovo(unqs[[i]], ..., verbose = verbose)
            }
            else if (method == "consensus") {
                bim <- isBimeraDenovoTable(unqs[[i]], ..., verbose = verbose)
            }
            else if (method == "per-sample") {
                bim <- t(apply(unqs[[i]], 1, function(x) isBimeraDenovo(x, 
                  ..., verbose = verbose)))
            }
            else {
                stop("Valid values for method: 'pooled', 'consensus', or 'per-sample'")
            }
            if (method %in% c("pooled", "consensus")) {
                outs[[i]] <- unqs[[i]][, !bim, drop = FALSE]
            }
            else if (method %in% c("per-sample")) {
                outs[[i]] <- unqs[[i]]
                outs[[i]][which(bim, arr.ind = TRUE)] <- 0
                cbim <- colSums(outs[[i]]) == 0
                outs[[i]] <- outs[[i]][, !cbim, drop = FALSE]
            }
            else {
                stop("Valid values for method: 'pooled', 'consensus', or 'per-sample'")
            }
        }
        else {
            stop("Unrecognized format: Requires named integer vector, dada-class, derep-class, sequence matrix, or a data.frame with $sequence and $abundance columns.")
        }
    }
    names(outs) <- names(unqs)
    if (length(outs) == 1) {
        return(outs[[1]])
    }
    return(outs)
}