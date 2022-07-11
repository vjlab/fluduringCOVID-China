# Tutorial 2 - The `baltic` tree object

This tutorial briefly covers the `tree`, `node` and `leaf` classes, and some of their most important attributes. Gytis' original Python2 `baltic` [readme](https://github.com/blab/baltic) is highly recommended to get an in-depth understanding on how the `bt.tree()` was set up.

![Image of baltic tree diagram](https://github.com/Don86/baltic3/blob/master/assets/baltic-tree-example.png)

* `tree` class: well, a tree. Retrieve all objects (internal nodes and leaves) as a list with `tre.Objects`
* `node` class: internal nodes of the tree.
    - `nodes`, annoyingly, have no `name` attribute. They can be given unique identifiers according to, say, enumerating by in-order traversal.[1]
    - A list of all nodes can be accessed via `tree.nodes`
    - Node x and y coordinates on the `matplotlib ax` can be retrieved with `node.height` and `node.y` respectively.
    - For chronograms, `node.absoluteTime` will have a computed date for internal nodes.
* `leaf` class: leaves of the tree.
    - The original sequence name used to compute the tree can be retrieved from `leaf.name`
    - Leaf x and y coordinates on the `matplotlib ax` can be retrieved with `leaf.height` and `leaf.y` respectively.
    - A list of all leaves can be accessed via `tree.leaves`
    - For chronograms, `leaf.absoluteTime` will have a computed date for leaves.


## Footnotes

[1] Even if they did have unique identifiers, non-visual subsetting of trees remains a problem (the main reason why I've wanted unique internal node identifiers has usually been to simply retrieve a subtree just from that internal node). But this can already be done in Figtree. If the tree gets too large to display in Figtree, it would similarly be too difficult to retrieve the internal nodes on a command line anyway. The most workable current solution is to donate to [icytree](https://icytree.org) for the developers to buy more server bandwidth, so that users can easily browse trees visually.
