#' Prioritise candidate variants by rarity and functional predictions
#' 
#' @param candidates candidate variants data.frame
#' @return candidate variants data.frame sorted with using very simple prioritisation scoring system
# #' @examples
# #' ***TODO***


# Prioritise candidate variants by rarity and predicted damage
# 
# Input:
#     candidates:        candidates data.frame (of type returned by read_candidate_variants or filter_candidate_variants functions)
#
# Returns:
#     prioritise_candidates data.frame

prioritise_candidates <- function(candidates)
{
    # *** VERY BASIC FUNCTION FOR PRIORITISING CANDIDATE VARIANTS -- UPDATE THIS IN FUTURE ***

    comphets <- endsWith(candidates$"inheritance model", "comp het")

    cand_comphets <- candidates[comphets, ]
    candidates <- candidates[!comphets, ]

    # Convert MAF to a -log10 score
    mean_MAF <- rowMeans(candidates[, c("MAF ExAC", "MAF ESP6500", "MAF 1000G")])
    MAF_score <- -log10(mean_MAF + 1e-8)

    # Convert damage predictions to a score, scaled so that it max damage == one order of magnitude rarer MAF with min damage
    SIFT_score <- as.numeric(gsub(".", 0, gsub("tolerated", -2, gsub("damaging", 2, candidates$SIFT)), fixed=TRUE))
    PP2_score <- as.numeric(gsub(".", 0, gsub("benign", -2, gsub("probably damaging", 2, gsub("possibly damaging", 1, candidates$Polyphen2))), fixed=TRUE))
    damage_score <- (SIFT_score + PP2_score) / 8

    RVIS_score <- (100 - ifelse(is.na(candidates$RVIS), 100, candidates$RVIS)) / 500

    prioritise_candidates <- candidates[order(MAF_score + damage_score + RVIS_score, decreasing=TRUE), ]
    prioritise_candidates <- rbind(prioritise_candidates, cand_comphets)

    return(prioritise_candidates)
}