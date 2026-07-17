# Shared helper: extract the numeric sample id from names like
# "RnDD_260713_12_RnDD-16S_n0" -> 12 (not the 260713 date field -- that one
# is followed by another digit, the real sample number is followed by a
# letter), so samples can be ordered/labelled 1, 2, 3 ... instead of
# alphabetically (1, 10, 11, 12 ...) or by full name.
# Falls back to alphabetical order (as-is) if the pattern doesn't match
# every sample name (e.g. non-RnDD naming).
sample_order_and_labels <- function(samples) {
    m <- regmatches(samples, regexpr("_([0-9]+)_[A-Za-z]", samples))
    nums <- suppressWarnings(as.integer(gsub("[_A-Za-z]", "", m)))
    if (length(m) == length(samples) && !any(is.na(nums))) {
        ord <- order(nums)
        list(order = ord, labels = as.character(nums)[ord])
    } else {
        ord <- order(samples)
        list(order = ord, labels = samples[ord])
    }
}
