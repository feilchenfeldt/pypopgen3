
import copy, json, os
import ete3
import dendropy

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt


def newick_to_node_name(nwk):
    """
    Create a formalized node name.
    """
    try:
        tree = HsTree(nwk)
    except ete3.parser.newick.NewickError:
        try:
            tree = HsTree(nwk + ';')
        except ete3.parser.newick.NewickError as e:
            raise e
    tree.sort_descendants()
    s = tree.write(format=9)[:-1]
    return s

class MassMigration(object):
    def __init__(self, source, destination, fraction, time):
        self.source = source
        self.destination = destination
        self.fraction = fraction
        self.time = time

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def get_time(self):
        return self.time

    @classmethod
    def from_dict(cls, dic, tree):
        source = tree.search_node_by_newick(dic['source_name'])
        destination = tree.search_node_by_newick(dic['destination_name'])
        mass_migration = cls(source, destination, dic['fraction'], dic['time'])
        return mass_migration

    def to_dict(self):
        return {'source_name':self.source.get_name(),
                'destination_name':self.destination.get_name(),
                 'fraction':self.fraction,
                 'time':self.get_time()}

class HsTree(ete3.Tree):
#    USE_NODE_DICT = True

    def __init__(tree, *args, **kwa):
        if args:
            newick = args[0]
        else:
            try:
                newick = kwa['newick']
            except KeyError:
                newick = None

        super(HsTree, tree).__init__(*args, **kwa)
        try:
            tree.mass_migrations = kwa['mass_migrations']
        except KeyError:
            tree.mass_migrations = []
            if newick is not None:
                if os.path.exists(newick+'.mass_migrations.json'):
                    mass_migration_ls = json.load(open(newick+'.mass_migrations.json'))
                    for ms in mass_migration_ls:
                        tree.mass_migrations.append(MassMigration.from_dict(ms, tree))



#        if HsTree.USE_NODE_DICT:
 #           tree._create_node_dict()

#    def _create_node_dict(tree):
#        max_dist = max([tree.get_distance(l) for l in tree.get_leaves()])
#        tree.node_dict = {}
#        for node in tree.traverse():
#            node.time  = max_dist - tree.get_distance(node)
#            node.name = node.get_name()
#            tree.node_dict.update({node.name: node})

#    def get_node(tree, node_name):
#        return tree.node_dict[node_name]

    @classmethod
    def from_dendropy(cls, dendro_tree):
        tree_str = dendro_tree.as_string(schema='newick').strip().split()[1]
        tree = cls(tree_str)
        #make sure no additional outer quotes in names with whitespace
        for t in tree.get_leaves():
                t.name = t.name.strip("' ")
        return tree

    def get_time(node):
        rt = node.get_tree_root()
        max_dist = rt.get_farthest_leaf()[1]
        return max_dist - rt.get_distance(node)

    def get_name(tree):
        node1 = copy.deepcopy(tree)
        node1.sort_descendants()
        s = super(HsTree, node1).write(format=9)[:-1]
        #s = node1.write(format=9)[:-1]
        return s



    def write(tree, **kwa):
        """
        ete3 write method with the addition
        of writing mass migrations into
        separate json.

        :return: (ete3 write output, mass_migrations)

        """

        output = super().write(**kwa)
        mass_migration_ls = []
        for mm in tree.mass_migrations:
            md = mm.to_dict()
            mass_migration_ls.append(md)

        try:
            mm_filename = kwa['outfile'] + '.mass_migrations.json'
            with open(mm_filename, 'w') as f:
                json.dump(mass_migration_ls, f)
        except KeyError:
            pass

        return output, mass_migration_ls


    def add_mass_migration(tree, source, destination, fraction, time):
#        if HsTree.USE_NODE_DICT:
#            source = tree.node_dict[source_name]
#            destination = tree.node_dict[destination_name]
#        else:
#            for node in tree.traverse():
#                if node.get_name() == source_name:
#                    source = node
#                if node.get_name() == destination_name:
#                    destintion = node
        mm = MassMigration(source, destination, fraction, time)
        if mm not in tree.mass_migrations:
            tree.mass_migrations.append(mm)

    def add_property_to_nodes(tree, property_name, property_node_dict):
        """
        Adds the attribute property_name to nodes with newick
        given as property_node_dict keys and property values as
        dictionary values.

        TODO: Update for case if node_dict is available.

        Example:
        print tree

                 /-A1
              /-|
           /-|   \-A2
          |  |
        --|   \-B
          |
           \-C

        property='ne'
        property_node_dict={'(A2,A1);':4, 'C;':1}

        This adds the property node.ne to the nodes:
           /-A1
        --|
           \-A2
        --C
        """
        dic = {}
        for k,v in property_node_dict.iteritems():
            dic[newick_to_node_name(k)] = v

        for node in tree.traverse():
            node_name = node.get_name()
            try:
                setattr(node, property_name, dic[node_name])
            except KeyError:
                pass


    def add_properties_to_nodes(tree, properties, properties_node_dict):
        """
        Adds the attributes in the list properties to nodes with newick
        given as property_node_dict keys and a dictionary of
        {property:value} as dictionary values.

        Example:
        print tree

                 /-A1
              /-|
           /-|   \-A2
          |  |
        --|   \-B
          |
           \-C

        properties=['ne', 'color']
        property_node_dict={'(A2,A1);': {'ne':4, 'color': 'black'},
                             'C;': {'color': 'green'}}

        """
        dic = {}
        for k, v in properties_node_dict.iteritems():
            dic[newick_to_node_name(k)] = v


        for node in tree.traverse():
            node_name = node.get_name()
            for prop in properties:
                try:
                    setattr(node, prop, dic[node_name][prop])
                except KeyError:
                    pass

    def get_nodes_at_time(tree, time):
        nodes = []
        for n in tree.traverse():
            ancestors = n.get_ancestors()
            if ancestors:
                if ancestors[0].get_time() > time > n.get_time():
                    nodes.append(n)
        return nodes


    def plot(tree, ax=None, style='orthogonal',
             orientation='left_to_right',
             ax_pos=None,
             origin=0,
             use_distances=True,
             leaves_to_present=False,
             internal_node_to_present=False,
             node_name_fun=None,
             node_name_format_fun=None,
             leaf_name_fun=None,
             leaf_name_format_fun=None,
             line_format_fun=None,
             migration_arrow_format_fun=None,
             tick_label_format_fun=None,
             debug=False):
        """
        Plot ete tree.
        """

        default_node_format_args = dict(xycoords='data', ha='center',
                                        xytext=(0, 1),
                                        textcoords='offset points',
                                        va='bottom',
                                        bbox=dict(boxstyle="round,pad=0.05", fc="w", alpha=0.5, lw=0),
                                        size=11)
        default_leaf_format_args = {'textcoords': 'offset points'}
        # , 'fontname':'monospace'
        default_line_format_args = {'color': 'k'}
        default_migration_arrow_format_args = dict(arrowstyle="->, head_length = 0.5, head_width = .5",
                                                   color='r', linestyle='solid', linewidth=2,
                                                   zorder=-1)

        if not use_distances:
            tree = HsTree(tree.write(format=9)[0])

        if internal_node_to_present:
            leaf_order = [n.get_name() for n in tree.traverse('preorder')]
        else:
            leaf_order = tree.get_leaf_names()

        if ax is None:
            fig = plt.figure(figsize=(12, len(leaf_order) * 0.3))
            ax = plt.gca()

        assert style in ['orthogonal', 'diagonal']
        assert orientation in ['left_to_right', 'bottom_to_top',
                               'top_to_bottom', 'right_to_left']

        if orientation == 'left_to_right':
            specific_default_leaf_format = {'va': 'center',
                                            'ha': 'left'}
            edgeline = ax.hlines
            orthogonalline = ax.vlines
            sign = 1
            order = 1
            ax_pos1 = 'bottom'


        elif orientation == 'bottom_to_top':
            specific_default_leaf_format = {'va': 'bottom',
                                            'ha': 'center',
                                            'rotation': 90}
            edgeline = ax.vlines
            orthogonalline = ax.hlines
            sign = 1
            order = -1
            ax_pos1 = 'right'


        elif orientation == 'top_to_bottom':
            specific_default_leaf_format = {'va': 'top',
                                            'ha': 'center',
                                            'rotation': 90}
            edgeline = ax.vlines
            orthogonalline = ax.hlines
            sign = -1
            order = -1
            ax_pos1 = 'right'


        elif orientation == 'right_to_left':
            specific_default_leaf_format = {'va': 'center',
                                            'ha': 'right'}
            edgeline = ax.hlines
            orthogonalline = ax.vlines
            sign = -1
            order = 1
            ax_pos1 = 'bottom'
        else:
            raise ValueError("Orientation must be left_to_right,"
                             "bottom_to_top, top_to_bottom or right_to_left.")

        if ax_pos is None:
            ax_pos = ax_pos1

        default_leaf_format_args.update(specific_default_leaf_format)

        default_leaf_format_args.update({'xytext': (sign * 5, 0)[::order]})

        # don't plot node names if no function given
        if node_name_fun is None:
            node_name_fun = lambda node: False
        if node_name_format_fun is None:
            node_name_format_fun = lambda node: {}
        # plot leaf.name as leaf name by default
        if leaf_name_fun is None:
            leaf_name_fun = lambda node: node.name
        if leaf_name_format_fun is None:
            leaf_name_format_fun = lambda node: {}
        if line_format_fun is None:
            line_format_fun = lambda node: {}
        if migration_arrow_format_fun is None:
            migration_arrow_format_fun = lambda node: {}
        if tick_label_format_fun is None:
            tick_label_format_fun = lambda x, p: format(-sign * int(x), ',')

        max_depth = tree.get_farthest_leaf()[1]

        max_label_width = 0
        leaf_annots = []

        for i, node in enumerate(tree.traverse('postorder')):
            time = node.get_time()
            if node.is_leaf():
                if leaves_to_present:
                    time = 0

                node.y = origin - leaf_order.index(node.name)

                leaf_name = leaf_name_fun(node)
                if leaf_name:
                    leaf_format_args = copy.deepcopy(default_leaf_format_args)
                    leaf_format_args.update(leaf_name_format_fun(node))
                    x = ax.annotate(leaf_name, xy=(-sign * time, sign * node.y)[::order],
                                    xycoords='data', **leaf_format_args)
                    leaf_annots.append(x)


            else:
                l = node.children[0]
                r = node.children[1]

                if not internal_node_to_present:
                    node.y = (l.y + r.y) / 2.
                else:
                    node.y = origin - leaf_order.index(node.get_name())

                # print(node.y)

                for c in (l, r):
                    line_format_args = copy.deepcopy(default_line_format_args)
                    line_format_args.update(line_format_fun(c))
                    if style == 'orthogonal':
                        if c.is_leaf() and leaves_to_present:
                            ctime = 0
                        else:
                            ctime = c.get_time()

                        edgeline(sign * (c.y), -sign * time, -sign * ctime, **line_format_args)
                        orthogonalline(-sign * time, *sorted([sign * c.y, sign * node.y]), **line_format_args)

                    elif style == 'diagonal':
                        ax.plot([-time, -c.get_time()], [node.y, c.y])

                    if not c.is_leaf():
                        node_name = node_name_fun(c)
                        if node_name:
                            node_format_args = copy.deepcopy(default_node_format_args)
                            node_format_args.update(node_name_format_fun(c))
                            ax.annotate(node_name, xy=((-time - c.get_time()) / 2., c.y),
                                        **node_format_args)

        for mm in tree.mass_migrations:
            # print "plotting migration one", mm.time, mm.source.get_name(), mm.destination.get_name()
            # ax.plot([-mm.time, -mm.time],sorted([mm.source.y, mm.destination.y]), color='r')
            # ax.arrow(-mm.time, mm.destination.y, 0 , mm.source.y - mm.destination.y,
            #                     length_includes_head=True, color='r', linestyle='dashed')
            migration_arrow_format_args = copy.deepcopy(default_migration_arrow_format_args)
            migration_arrow_format_args.update(migration_arrow_format_fun(c))
            ax.annotate("", xytext=(-sign * mm.time, sign * mm.destination.y)[::order],
                        xy=(-sign * mm.time, sign * mm.source.y)[::order],
                        arrowprops=migration_arrow_format_args)

            ax.annotate("{}%".format(int(round(mm.fraction * 100))),
                        xy=(-sign * mm.time, sign * (mm.destination.y + mm.source.y) / 2.)[::order],
                        xytext=(-sign * 5, 0)[::order],  # ha='right',va='center',
                        bbox=dict(boxstyle="round,pad=0.1", fc="w", alpha=0.5, lw=0),
                        textcoords='offset points', color='r')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        if ax_pos:
            ax.spines[ax_pos].set_visible(True)
        else:
            if not debug:
                ax.set_yticks([])
                ax.set_xticks([])
                plt.tick_params(
                    # axis='x',  # changes apply to the x-axis
                    which='both',  # both major and minor ticks are affected
                    left=False,  # ticks along the bottom edge are off
                    right=False,  # ticks along the top edge are off
                    bottom=False,
                    top=False,
                    labelbottom=False)  # labels along the bottom edge are off

        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()

        if ax_pos in ['top', 'bottom']:
            ax.set_ylim([ymin - (ymax - ymin) * 0.05, ymax + (ymax - ymin) * 0.01])
            # ax.xaxis.tick_bottom()
            if not debug:
                ax.set_yticks([])
            ax.get_xaxis().set_major_formatter(
                mpl.ticker.FuncFormatter(tick_label_format_fun))

        elif ax_pos in ['left', 'right']:

            ax.set_xlim([ymin - (ymax - ymin) * 0.05, ymax + (ymax - ymin) * 0.01])
            if not debug:
                ax.yaxis.tick_right()
                ax.set_xticks([])
            ax.get_yaxis().set_major_formatter(
                mpl.ticker.FuncFormatter(tick_label_format_fun))

        fig = plt.gcf()
        renderer = fig.canvas.get_renderer()
        inv = ax.transData.inverted()

        labeledges = []

        for x in leaf_annots:
            bb = x.get_window_extent(renderer=renderer)
            ((x0, y0), (x1, y1)) = inv.transform(bb)
            labeledges.append(((x0, y0), (x1, y1)))

        labeledges = np.array(labeledges)
        # print(np.sort(labeledges[:, int((sign + 1) / 2), int((order - 1) / -2)])[::-order])

        edgecoord = np.sort(labeledges[:, int((sign + 1) / 2), int((order - 1) / -2)])[::-sign][0]
        # print(edgecoord)

        if internal_node_to_present:
            for node in tree.traverse():
                if not node.is_leaf():
                    edgeline(sign * node.y, -sign * node.get_time(),
                             edgecoord,
                             linestyle='dotted')  # dotted

        # edgecoord = edgecoord * 1.1

        # if orientation == 'left_to_right':
        #     ec = max(xmax, edgecoord ) + (xmax-xmin)*0.05
        #     #ax.set_xlim(xmin, ec)
        #     ax.annotate('Y',xy=(edgecoord, 0),
        #                ha='left')
        #     plt.scatter(edgecoord, 0)
        #     print(labeledges,edgecoord)
        # elif orientation == 'right_to_left':
        #     ec = min(xmin, edgecoord) - (xmax-xmin)*0.05
        #     ax.set_xlim(ec,xmax)
        # elif orientation == 'top_to_bottom':
        #     ax.set_ylim(min(ymin, edgecoord), ymax)
        # elif orientation == 'bottom_to_top':
        #     ax.set_ylim(ymin, max(ymax, edgecoord ))

        # ymin, ymax = ax.get_ylim()
        # ax.set_ylim([ymin-(ymax-ymin)*0.05,ymax+(ymax-ymin)*0.01])
        # ax.xaxis.tick_bottom()
        # ax.get_xaxis().set_major_formatter(
        #        mpl.ticker.FuncFormatter(xtick_label_format_fun))
        # ax.set_xlim([1,10])

        return ax, leaf_annots

    def search_node_by_newick(tree, newick):
        for n in tree.traverse():
            if n.get_name() == newick:
                return n
        raise Exception('Node not found.')

    def set_leaf_order(tree, order, check_consistent=True):
        """
        Changes the tree so that the leaves
        conform to order.
        The order must be consistent with
        the branching structure.

        Parameters:
        tree ... ete3 tree object
        order ... list of leaf names

        Returns:
        None (tree changed in place)
        """
        order = list(order)
        for i, node in enumerate(tree.traverse('postorder')):
            if not node.is_leaf():
                l = node.children[0]
                r = node.children[1]
                lnames = l.get_leaf_names()
                rnames = r.get_leaf_names()
                if order.index(lnames[0]) > order.index(rnames[0]):
                    node.swap_children()
        if check_consistent:
            if tree.get_leaf_names() != order:
                raise Exception("The provided order is not conistent with tree: \nResulting order: {}\nInput order:{}".format(
                    tree.get_leaf_names(), order))


    def set_outgroup(tree, outgroup, end_at_present=True):
        """
        Set root of tree.

        :param outgroup: name of the outgroup leaf
        :param end_at_present: root at a point of the branch so that outgroup ends at time 0.
        :return: None
        """
        super().set_outgroup(outgroup)
        if end_at_present:
            outgroup_node = tree.search_nodes(name=outgroup)[0]
            ingroup_root = [n for n in tree.get_children() if n is not outgroup_node][0]
            time = outgroup_node.get_time()
            outgroup_node.dist = outgroup_node.dist + time / 2
            ingroup_root.dist = ingroup_root.dist - time / 2
            assert ingroup_root.dist > 0, \
                    "Outgroup branch too short to lead to present. Use end_at_present=False."

    def reverse(tree):
        tree.set_leaf_order(tree, tree.get_leaf_names()[::-1])


def get_random_geneflow_species_tree(n_species,
                                     n_geneflow_events,
                                     max_depth=1e6,
                                     effective_population_sizes=5e4,
                                     sample_sizes=2,
                                     all_branches_end_at_present=True,
                                     add_outgroup=True,
                                     outgroup_dist_factor=2,
                                     geneflow_times_rand_fun=np.random.uniform,
                                     geneflow_strength_rand_fun=lambda: np.random.beta(2, 5) * 0.3,
                                     random_tree_process='birth_death'):
    """
    This function generates a random HsTree tree of n_species,
    with n_geneflow_events random, unidirectional gene flow (MassMigration)
    events distributed along the tree.
    The resulting species tree can be used as an input for msprime coalescent
    simulations, by using the pypopgen3 function simulate.simulate_from_tree.

    The tree is randomly generated using ete3's Tree.populate function.
    Gene-flow times are drawn randomly, distribution can be provided (default:uniform).
    Gene-flow partners are drawn with equal probablitites form all branches
    that exist at the gene-flow time (except outroup).
    Gene-flow strengths are drawn from given distribution (default: np.random.beta(2, 5) * 0.3)


    :param n_species: Number of (ingroup) species.
    :param n_geneflow_events: Number of gene flow (MassMigration) events.
    :param max_depth: depth (in generations) of the ancestral node of ingroup species.
    :param effective_population_sizes: Effective population sizes of ingroup species. (default: 1e6)
                                      Can be integer or list of length n_species. (default: 5e4)
    :param sample_sizes: Sample sizes (2n) of ingroup species (for msprime).
                        Can be integer or list of length n_species. (default: 2)
                        Should be in multiples of 2 because we assume diploid organisms.
    :param all_branches_end_at_present: If true, all branches end at the same time. (default: True)
    :param add_outgroup: If true, add an outgroup with sample size 2n=2 and Ne=1.
                        Hence, the outgroup is not variable. (default: True)
    :param outgroup_dist_factor: Factor by which to multiply max_depth,
                                to get the divergence of the outgroup. (default: 2)
    :param geneflow_times_rand_fun: Random function to determine gene flow times.
                                    (default: np.random.uniform)
    :param geneflow_strength_rand_fun: Random function to determine gene flow strength.
                                    (default: np.random.beta(2, 5) * 0.3)
    :param random_tree_process: Process that creates the random tree. Possible values are
                                'birth_death' ... yule process tree created with dendropy
                                                  with death rate 0 and birth rate to match
                                                  max_depth and number of species
                                                  This corresponds to the expectation for an actual
                                                  species tree
                                'ete3' ... tree created with ete3 and random branch lenght
                                            This looks more like a coalescent tree with comparatively
                                            long terminal branches.
    :return:
        HsTree object (derived from ete3 tree)

    """
    tree_processes = ['birth_death', 'ete3']

    assert random_tree_process in tree_processes

    try:
        effective_population_sizes = list(effective_population_sizes)
        assert len(effective_population_sizes) == n_species
    except TypeError:
        effective_population_sizes = [effective_population_sizes] * n_species

    try:
        sample_sizes = list(sample_sizes)
        assert len(sample_sizes) == n_species
    except TypeError:
        sample_sizes = [sample_sizes] * n_species

    if random_tree_process == 'ete3':

        t = HsTree()
        t.populate(n_species, random_branches=True)

        # # Make as short as possible unique names
        # # Remove leading uninformative aaaa...
        # name_arr = np.array([list(n) for n in t.get_leaf_names()])
        # for i in range(name_arr.shape[1]):
        #
        #     if len(set(name_arr[:, i])) > 1:
        #         name_informative_index = i
        #         break
        # for l, n in zip(t.iter_leaves(), name_arr[:, name_informative_index:]):
        #     l.name = "".join(n)


    elif random_tree_process == 'birth_death':
        t0 = dendropy.simulate.treesim.birth_death_tree(birth_rate=1, death_rate=0,
                                                       num_extant_tips=n_species, gsa_ntax=n_species+1)
        t = HsTree.from_dendropy(t0)
    else:
        raise ValueError(f"random_tree_process must be in {tree_processes}")

    digits = int(np.ceil(np.log10(n_species)))
    for i, n in enumerate(t.iter_leaves()):
        n.name = 'S' + "{{:0{}d}}".format(digits).format(i)


    depth = t.get_time()

    if all_branches_end_at_present:
        for l in t.iter_leaves():
            root_to_node = t.get_distance(l)
            l.dist = l.dist + depth - root_to_node

    for ne, n_samples, node in zip(effective_population_sizes,
                                   sample_sizes, t.get_leaves()):
        node.ne = ne
        node.n_samples = n_samples

    for n in t.traverse(strategy='postorder'):
        n.dist = n.dist * max_depth / depth
        if not n.is_leaf():
            nes = [sn.ne for sn in n.get_children()]
            n.ne = sum(nes) * 1. / len(nes)

    depth = t.get_time()

    if add_outgroup:
        t2 = HsTree()

        t2.add_child(t)
        t.dist = max_depth * (outgroup_dist_factor - 1)
        outgroup_node = t.add_sister(name='outgroup', dist=max_depth * outgroup_dist_factor)
        outgroup_node.ne = 1
        outgroup_node.n_samples = 2

    else:
        t2 = t

    for i in range(n_geneflow_events):
        time = geneflow_times_rand_fun() * depth
        nodes_at_time = t.get_nodes_at_time(time)
        source, destination = np.random.choice(nodes_at_time, 2, replace=False)
        fraction = geneflow_strength_rand_fun()
        t.add_mass_migration(source, destination, fraction, time)

    t2.mass_migrations = t.mass_migrations

    return t2


HsTree.write.__doc__ =  HsTree.write.__doc__ + ete3.Tree.write.__doc__


def get_fmin_tree(f_df, tree):

    """
    """

    f = f_df[f_df['F4ratio']>=0].reset_index()

    t = copy.deepcopy(tree)

    i=0
    for node in  t.traverse():
        if node.children:
            l = node.children[0]
            r = node.children[1]
            lleaves = l.get_leaf_names()
            rleaves = r.get_leaf_names()

            node_fl = f[f['h2'].isin(lleaves)&f['h1'].isin(rleaves)]
            node_fr = f[f['h2'].isin(rleaves)&f['h1'].isin(lleaves)]



            for side, node_f, sister_f in [(0,node_fl, node_fr),(1,node_fr, node_fl)]:
                if len(node_f) or len(sister_f):

                    sister_f0 = sister_f.rename(columns={'h1':'h2','h2':'h1'})
                    sister_f0['F4ratio'] = 0
                    sister_f0['Z'] = 0
                    nf = pd.concat([node_f, sister_f0])

                    #node_f.sort_values('|f|', ascending=False)
                    #only take h3 with maximum mean '|f|' on this branch
                    #h3 = node_f.groupby('h3').mean().sort_values('|f|', ascending=False).iloc[0].name
                    #node_f1 = node_f[node_f['h3']==h3]
                    child = node.get_children()[side]

                    #child.add_feature('rscore', summary(node_f1['|f|']))
                    #child.add_feature('h3', h3)
                    child.add_feature('branch_f', nf.groupby(['h2','h3']).min().reset_index())

    return t

def get_node_name(node):
    if node.is_leaf():
        return node.name
    else:
        return str(",".join(["".join([n[0] for n in c.get_leaf_names()]) for c in node.get_children()]))

def try_get_f(node, taxa, statistic='F4ratio', cp_summary=np.nanmean):
    if hasattr(node, 'branch_f'):
        h3groups = node.branch_f.groupby('h3')
        h3_summary = h3groups.apply(lambda df: cp_summary(df[statistic].values))
        return h3_summary.ix[taxa]
    else:
        return pd.Series({t:np.nan for t in taxa})

def get_branch_mat(rscore_tree, statistic='F4ratio',cp_summary=np.nanmean):
    """
    Tree without outgroup.
    """
    taxa = rscore_tree.get_leaf_names()
    branch_mat_df = pd.DataFrame()
    for node in rscore_tree.iter_descendants('preorder'):
        node_name = get_node_name(node)
        while node_name in branch_mat_df.columns:
            node_name = node_name + '0'
        branch_mat_df.loc[:,node_name] = try_get_f(node, taxa, statistic=statistic, cp_summary=cp_summary)#cp_summary=np.nanmax
    branch_mat = branch_mat_df.T.loc[:,taxa]
    branch_mat = branch_mat.iloc[::-1,:]
    return branch_mat



def align_fbranch_with_tree(fbranch, tree, outgroup, ladderize=False):
    tree_no = copy.deepcopy(tree)
    # tree_no.ladderize()
    # remove outgroup
    tree_no.prune([n for n in tree_no.get_leaf_names() if n != outgroup])
    #
    if ladderize:
        tree_no.ladderize()

    fb = fbranch.copy()
    fb.index = fb.index.droplevel(0)
    try:
        fb = fb.drop(outgroup, axis=1)
    except KeyError:
        pass
    fb.index = [tuple(sorted(i.split(','))) for i in fb.index]
    row_order = []
    col_order = tree_no.get_leaf_names()
    for n in tree_no.iter_descendants(strategy='preorder'):
        row_order.append(tuple(sorted(n.get_leaf_names())))

    assert fb.shape == (len(row_order), len(col_order)), \
        "{} != ({},{}))".format(fb.shape, len(row_order), len(col_order))

    assert set(fb.index.values) == set(row_order), \
        "Samples in input tree and fbranch matrix not consistent: {} != {}".format(fb.index.values, row_order)

    assert set(fb.columns.values) == set(col_order), \
        "Samples in input tree and fbranch matrix not consistent: {} != {}".format(fb.columns.values, col_order)

    # order fb in a way that is conistent with the tree
    fb = fb.loc[row_order, col_order].iloc[::-1]  # ,::-1

    return fb, tree_no



def plot_fbranch(fbranch, tree_no_outgroup, leaves_to_present=True,
                 use_distances=False,
                 debug=False, tree_label_size=14, max_color_cutoff=None):
    # print("1706")

    n_rows, n_cols = fbranch.shape

    if max_color_cutoff is None:
        max_color_cutoff = fbranch.max().max()

    #plt.rcParams['font.size'] = 12
    plt.rcParams['font.size'] = tree_label_size


    # depth = tree_sd2_1_ete.get_farthest_leaf(topology_only=True)[1] +2

    # visited = []

    # fig = plt.figure(figsize=(18,20))

    fig = plt.figure(figsize=(n_cols, n_rows * 0.5))

    toptree_ax = plt.subplot2grid((5, 7), (0, 2), rowspan=1, colspan=4)

    tree_no_outgroup.plot(origin=-0.5, ax=toptree_ax,
                          orientation='top_to_bottom',
                          use_distances=use_distances,
                          leaves_to_present=leaves_to_present)

    if not debug:
        toptree_ax.axis('off')

    lefttree_ax = plt.subplot2grid((5, 7), (1, 0), rowspan=4, colspan=2)

    tree_no_outgroup.plot(origin=fbranch.shape[0] + 0.5, ax=lefttree_ax, orientation='left_to_right',
                          use_distances=use_distances,
                          leaves_to_present=leaves_to_present,
                          internal_node_to_present=True)
    if not debug:
        lefttree_ax.axis('off')
    # ax = plot_node_tree(tree_no_outgroup, ax=ax, x0=0,y0=len(branch_mat)-0.5,em=0.5,fontsize=12)

    # ax.set_xlabel('excess allele sharing with')
    # toptree_ax.xaxis.set_label_position('top')

    # ax.set_xlim([-1,10])

    fbranch_ax = plt.subplot2grid((5, 7), (1, 2), rowspan=4, colspan=4, sharey=lefttree_ax, sharex=toptree_ax)

    branch_mat0 = fbranch.copy()
    branch_mat0_masked = np.ma.array(branch_mat0, mask=np.isnan(branch_mat0))

    fmax = branch_mat0.max().max()
    fmin = branch_mat0.min().min()
    colors = np.concatenate([[[1, 1, 1, 1]], plt.cm.Reds(np.linspace(0., 1 * fmax / max_color_cutoff, 256))])
    mymap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)

    plt.pcolormesh(branch_mat0_masked, cmap=mymap, rasterized=True)  # ,cmap=jet

    # plt.scatter(zs.index.droplevel(0).values+0.5,zs.index.droplevel(1)+0.43,marker='*',s=10)

    fbranch_ax.set_xticks(np.arange(0.5, fbranch.shape[1]))
    # fbranch_ax.set_xticklabels(fbranch.columns,rotation=90)
    fbranch_ax.set_xticklabels([])
    fbranch_ax.set_yticks(np.arange(0.5, fbranch.shape[0]))
    fbranch_ax.set_yticklabels([])
    fbranch_ax.set_facecolor((0.85, 0.85, 0.85))
    # plt.tight_layout()
    fbranch_ax.set_ylim(0, len(fbranch))
    # plt.subplots_adjust(wspace=-0.1)
    # plt.subplots_adjust(hspace=0.05)
    fbranch_ax.xaxis.tick_top()
    fbranch_ax.set_xlim([0, fbranch.shape[1]])

    for b in range(fbranch.shape[0]):
        l = fbranch_ax.axhline(y=b + 0.04, xmin=0, xmax=1, linewidth=1, color='grey', alpha=0.5)

    for b in range(fbranch.shape[1]):
        l = fbranch_ax.axvline(x=b - 0.02, ymin=0, ymax=1, linewidth=1, color='grey', alpha=0.5)

    plt.tight_layout()

    cb_ax = fig.add_axes([fbranch_ax.get_position().xmax + 0.02, fbranch_ax.get_position().ymin,
                          0.03, fbranch_ax.get_position().ymax - fbranch_ax.get_position().ymin])

    # plt.subplot2grid((5,7), (1,6), rowspan=4,colspan=1)
    mappable = mpl.cm.ScalarMappable(cmap=mymap)
    mappable.set_array([fmin, fmax])
    cbar = plt.colorbar(mappable, cax=cb_ax, label='$f_b$')
    cbar.set_label('$f_b$', size=20 * (1 + np.log(n_cols / 8.)))
    cbar.solids.set_rasterized(True)
    cbar.ax.tick_params(labelsize=14 * (1 + np.log(n_cols / 8.)))

    # plt.axis('tight')

    # plt.subplots_adjust(wspace=0.05)

    # plt.tight_layout()

    return toptree_ax, lefttree_ax, fbranch_ax, cb_ax
