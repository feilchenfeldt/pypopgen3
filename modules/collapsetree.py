from collections import Counter
import copy
import numpy as np
from treetools import HsTree as HsTree

def get_collapsed_tree(ind_tree, sample_to_species):
    """
    Collapses leafs in a tree for all monophyletic groups
    of samples from the same species. If species are not
    monophyletic there will be mutliple leafs per species.

    :param ind_tree:
    :param sample_to_species: dict-like were key is sample id and value is species
    :return:
    """
    summary_tree = copy.deepcopy(ind_tree)
    sts = sample_to_species
    #already_collapsed = []
    to_collapse = []


    tree_samples = summary_tree.get_leaf_names()
    sample_groups = {}
    group_dists = {}

    for n, nn in zip(ind_tree.iter_leaves(),summary_tree.iter_leaves()):
        name = n.name
        if name in to_collapse:
            continue
        ids = [name]
        old_a = n
        for a in n.get_ancestors():
            new_ids = a.get_leaf_names()
            try:
                species_names = set([sts[nm] for nm in new_ids])
            except Exception as e:
                print(new_ids)
                raise e
            species_number = len(species_names)
            if species_number <= 1:
                ids = new_ids
                old_a = a

            else:
                species = sts[name]
                i = 0
                while species + '.' + str(i) in sample_groups.keys():
                    i += 1
                group_name = species + '.' + str(i)
                sample_groups.update({group_name:ids})
                #mean distance of samples to common ancestor
                if len(ids) > 1:
                    mean_dist = np.mean([nx.get_distance(old_a) for nx in old_a.get_leaves()])
                else:
                    mean_dist = n.dist
                to_collapse = to_collapse + [id for id in ids if id != name ]
                group_dists.update({group_name:mean_dist})
                #nn.dist = mean_dist
                nn.name = group_name
                break

            #print()
            #break
        #break
    summary_tree.prune(sample_groups.keys())
    for l in summary_tree.iter_leaves():
        l.dist =  group_dists[l.name]

    return summary_tree, sample_groups


class CollapsedTree(HsTree):
    def __init__(self, ind_tree, sample_to_species):
        collapsed_tree, sample_groups  = get_collapsed_tree(ind_tree, sample_to_species)
        self.tree = collapsed_tree
        self.groups = sample_groups
        self.sample_tree = ind_tree
        self.sample_to_species = sample_to_species