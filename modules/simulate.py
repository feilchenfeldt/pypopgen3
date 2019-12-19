"""
Tools to simulate data using
msprime.
"""
import subprocess
import msprime
import numpy as np


def msprime_input_from_split_tree(stree, diploid=True):
    """
    Get msprime input.
    Only works if nodes are annotated with ne and
    n_samples. E.g. if tree was produced with
    pypopgen.modules.splittree.
    """
    leaf_order = stree.get_leaf_names()
    name_to_id = {n: i for i, n in enumerate(leaf_order)}
    id_to_name = {v: k for k, v in name_to_id.items()}

    population_configurations = [0] * len(leaf_order)
    samples = []
    demographic_events = []
    event_times = []
    pop_to_samples = {}

    for i, node in enumerate(stree.traverse('postorder')):
        try:
            node.ne = float(node.ne)
        except AttributeError:
            pass
        try:
            node.n_samples = int(node.n_samples)
        except AttributeError:
            pass

        time = round(node.get_time(), 2)
        if node.is_leaf():
            if diploid:
                assert int(node.n_samples) % 2 == 0, \
                    ("If diploid=True, n_samples must be a multiple of 2, "
                     " i.e., the haploid number of samples. But for node {} "
                     " n_samples=={}".format(node.name, node.n_samples))
            node.id = name_to_id[node.name]
            population_configurations[node.id] = msprime.PopulationConfiguration(initial_size=node.ne)
            pop_to_samples[node.id] = []
            for j in range(int(node.n_samples)):
                samples.append(msprime.Sample(population=node.id, time=time))
            for k in range(int(int(node.n_samples) / 2) if diploid else int(node.n_samples)):
                pop_to_samples[node.id].append(node.name + '.' + str(k))
        else:
            l = node.children[0]
            r = node.children[1]
            node.id = l.id

            # Population merging
            event = msprime.MassMigration(
                time=time, source=r.id, destination=node.id,
                proportion=1.0)
            demographic_events.append(event)
            event_times.append(time)

            # Population size changes to mean of populations below
            if not node.is_root() and node.ne != l.ne:
                ne_change = msprime.PopulationParametersChange(
                    time=time, initial_size=node.ne, population_id=node.id)
                demographic_events.append(ne_change)
                event_times.append(time)

    for mm in stree.mass_migrations:
        introgression_event = msprime.MassMigration(
            time=mm.time,
            source=mm.source.id,
            destination=mm.destination.id,
            proportion=mm.fraction)
        demographic_events.append(introgression_event)
        event_times.append(mm.time)

    sample_names = [s for v in pop_to_samples.values() for s in v]
    sample_to_pop = {s: id_to_name[k] for k, v in pop_to_samples.items() for s in v}
    sorted_events = list(np.array(demographic_events)[np.argsort(event_times)])

    return samples, population_configurations, sorted_events, sample_to_pop, sample_names, id_to_name


def simulate(samples, population_configurations, demographic_events,
             mutation_rate, recombination_rate, genomic_length):
    tree_sequence = msprime.simulate(samples=samples,
                                     recombination_rate=recombination_rate,
                                     mutation_rate=mutation_rate,
                                     length=genomic_length,
                                     population_configurations=population_configurations,
                                     demographic_events=demographic_events)
    return tree_sequence


def simulate_from_tree(stree, mutation_rate, recombination_rate,
                       genomic_length, filename=None, chrom_id=1, save_vcf=True, diploid=True):
    if save_vcf:
        assert filename is not None, "provide filename if save_vcf is True"
    else:
        assert chrom_id == 1, "setting chrom_id is only supported if saved as vcf"
    samples, population_configurations, sorted_events, sample_to_pop, sample_names, id_to_name \
        = msprime_input_from_split_tree(stree, diploid=diploid)
    if not save_vcf:
        tree_sequence = simulate(samples, population_configurations, sorted_events,
                                 mutation_rate, recombination_rate, genomic_length)
        return id_to_name, sample_to_pop, sample_names, tree_sequence

    else:
        simulate_to_vcf(filename, samples, population_configurations, sorted_events,
                        mutation_rate, recombination_rate, genomic_length, sample_names=sample_names,
                        chrom_id=chrom_id, diploid=diploid)
        return id_to_name, sample_to_pop, sample_names, None


def simulate_to_vcf(filename, samples, population_configurations, demographic_events,
                    mutation_rate, recombination_rate, genomic_length, sample_names=None, chrom_id=1, diploid=True):
    tree_sequence = simulate(samples, population_configurations, demographic_events,
                             mutation_rate, recombination_rate, genomic_length)
    tree_sequence_to_vcf(filename, tree_sequence, sample_names, chrom_id, diploid)


def tree_sequence_to_vcf(filename, tree_sequence, sample_names=None, chrom_id=1, diploid=True):
    format_commands = ""
    if sample_names is not None:
        format_commands += ("awk ' {{if ($0 ~ \"^#CHROM\") {{$10=\"{}\"; "
                            "for (i=1; i<=10; i++) printf $i \"\t\"; print \"\"}} else print $0}}' | ").format(
            "\t".join(sample_names))
    #if chrom_id != 1:
    #    format_commands += "awk 'BEGIN{{OFS=\"\t\"}} {{if ($1 == 1) $1={};print $0}}' | ".format(chrom_id)
    format_commands += "bgzip -c"

    bgzip_stream = subprocess.Popen([format_commands],
                                    stdin=subprocess.PIPE,
                                    stdout=open(filename, 'w'),
                                    stderr=subprocess.PIPE, shell=True,
                                    encoding='utf8')
    tree_sequence.write_vcf(bgzip_stream.stdin, 2 if diploid else 1, contig_id=chrom_id)
    bgzip_stream.communicate()
    tabix_p = subprocess.Popen(['tabix', '-f', '-p', 'vcf', filename])
    tabix_p.wait()