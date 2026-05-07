from skbio import TreeNode

## Copy from https://github.com/scikit-bio/scikit-bio/blob/0.6.2/skbio/tree/_tree.py#L1371:
def root_at_midpoint(reset=False, branch_attrs=["name"], root_name="root") -> TreeNode:
    """Reroot the tree at the midpoint of the two tips farthest apart."""
        
    if reset:
        tree.unroot()

    max_dist, tips = tree.get_max_distance()
    half_max_dist = max_dist / 2.0

    if max_dist == 0.0:
        return tree

    tip1 = tree.find(tips[0])
    tip2 = tree.find(tips[1])
    lca = tree.lowest_common_ancestor([tip1, tip2])

    if tip1.accumulate_to_ancestor(lca) > half_max_dist:
        climb_node = tip1
    else:
        climb_node = tip2

    dist_climbed = 0.0
    while dist_climbed + climb_node.length < half_max_dist:
        dist_climbed += climb_node.length
        climb_node = climb_node.parent

    # case 1: midpoint is at the climb node's parent
    # make the parent node as the new root
    if dist_climbed + climb_node.length == half_max_dist:
        new_root = climb_node.parent

    # case 2: midpoint is on the climb node's branch to its parent
    # insert a new root node into the branch
    else:
        new_root = tree.__class__()
        climb_node.insert(new_root, half_max_dist - dist_climbed)
        # TODO: Here, `branch_attrs` should be added to `insert`. However, this
        # will cause a backward-incompatible behavior. This change will be made
        # in version 0.7.0, along with the removal of `name` from the default of
        # `branch_attrs`.

    branch_attrs = set(branch_attrs)
    branch_attrs.update(["length", "support"])
    return new_root.unrooted_copy(branch_attrs=branch_attrs, root_name=root_name)

if __name__ == "__main__":
    
    # Original
    tree = TreeNode.read(["((a:1,b:1)c:2,(d:3,e:4)f:5,g:1)h;"])
    print(tree.ascii_art())

    # Find tree root
    rooted_tree  = root_at_midpoint(tree)
    print(rooted_tree.ascii_art())


