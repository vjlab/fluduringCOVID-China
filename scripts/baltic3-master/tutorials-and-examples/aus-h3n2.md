This analysis is powered by [`treetime`](https://github.com/neherlab/treetime).

The figure is a *time-stamped phylogenetic tree*, or *chronogram*. The tips of the tree at the right side of the figure are the sequence samples, sampled at various points in time (roughly over the last 10 years). The internal nodes of the tree are hypothetical ancestors of these tips (and *their* respective hypothetical ancestors), so that we have some notion of the "family tree" of virus sequences that we have.

At about 2011, a rather disturbing thing started to happen: more lineages of the A/H3N2 virus started circulating at once. Before 2011, A/H3N2 viruses generally only had 1 "main" lineage. After that,



This next figure is a different kind of phylogenetic tree. This time, the branches aren't related to time, but genetic distance: a longer branch length indicates a greater genetic distance (whereas the previous figure placed the nodes in the tree according to their occurrence in *time*).

The `treetime` program has kindly computed for us the maximmum-likelihood estimate (translation: "best guess") of which state each ancestor originated from, given input information about the originating state of each virus sequence I gave it. The most obvious pattern here is that there is no pattern; viruses appear to be able to jump wildly from any state to any state. This is a bit of a puzzle, since Australia is a pretty damn big country; the flu virus may be able to hang around in the air a little bit, but it's a bit of a stretch to imagine that it can be blown on the wind from Sydney to Perth in less than a month. Other suggestions include viruses being carried by sick people by air or land travel.

## About the Data

## Technical Notes

[2] 
