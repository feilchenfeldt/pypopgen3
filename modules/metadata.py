def make_change(mt, sample, column, new_value, note=None):
    mt.loc[sample, "change_history"] = column + ':'  + str(mt.loc[sample, column]) + ';'
    mt.loc[sample, column] = new_value
    if note is not None:
        mt.loc[sample, 'notes'] = mt.loc[sample, 'notes'] + ';' + note
    if column in ['genus', 'species']:
        mt['full_name'] = mt['genus'] + '_' + mt['species']


def get_pub_species_name(s):
    return ("cf. " if s['taxonomic_uncertainty'] == 'cf.' else '') + ("sp. '" if s['undescribed'] else '') + s[
        'species'].replace('-', ' ') + ("'" if s['undescribed'] else '')

def update_names(mt):
    mt['full_name'] = mt['genus'] + '_' + mt['species']
    quotes = mt['undescribed'].apply(lambda s: "'" if s == 1 else "$")
    genus = mt['genus'].apply(lambda s: s[0])
    species = mt['species'].apply(lambda s: s.replace('-', ' '))
    mt['plot_name'] = "$" + genus + "$. " + ("cf. " if mt['taxonomic_uncertainty'] == 'cf.' else '')  + quotes + species + quotes






