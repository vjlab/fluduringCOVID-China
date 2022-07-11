# Tutorial 1 - Input Tree File Format
`Baltic` was originally written to accept BEAST MCC trees in mind, so any `.tre` file or `.nex` file needs to be wrangled into a `nexus` format. Most tree computation software like `raxml`, `fasttree` or `treesub` will only spit out the newick string. There are 2 blocks which need for be formatted: the taxa block, and the tree block.

Say, we have a newick string like this:

```
(A_2010:0.1,B_2017:0.2,(C_2015:0.3,D_2012:0.4)[&branch="B1"]:0.5);
```

We'll need to add some stuff around it, so that `baltic` knows where to read in tipnames, labels, and so on (1).

### "Taxa" Block

This is a list of the taxa/virus/tip names, `Taxa1`, `Taxa2`...etc.

```
Begin taxa;
    Dimensions ntax=300;
        Taxlabels
        Taxa1
        Taxa2
        ...
        Taxa300
;
End;
```

(Previously, the `austechia_read_tree()` function required tipnames to have dates appended to the end; that requirement has been removed.)

### "Trees" Block
This is the block containing the actual tree itself, in Newick format.

```
begin trees;
    tree TREE1 = [&R] <start of newick format...>
```

Altogether, the following `Nexus` format is a valid format. Note that the node labels, `[&branch="something"]`, are optional. These node/leaf proerties are usually some property pertaining to the software output that generated the tree, e.g. amino acid mutations, bootstrap values, etc, and can be accessed via the `node.traits` or `leaf.traits` attributes. If this is, for instance, a RAXML tree with bootstrap labels at the nodes, then it would be something like `[&label=93]`, or something like that. The following is a full example of a formatted input `.tre` file:

```
#NEXUS
begin taxa;
    dimensions ntax=4;
	Taxlabels
	A
	B
	C
	D
;
end;

begin trees;
    tree TREE1 = [&R] (A:0.1,B:0.2,(C:0.3,D:0.4)[&branch="B1"]:0.5)[&branch="B2"];
end;
```

Produces a tree in Baltic tree like:

![Image of Wikitree](https://github.com/Don86/baltic3/blob/master/assets/my_simple_tre.png)
