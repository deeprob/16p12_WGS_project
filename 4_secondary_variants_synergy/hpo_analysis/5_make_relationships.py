#!/bin/python
import pandas as pd

# Make files with HPO relationships
obo = pd.read_csv('Analysis_files/3_hpo_obo_reformat.csv')

output = obo[['id', 'alt_id', 'is_a']].copy()
# Add "direct" term relationships
# A direct relationship looks like this:
# If the hierarchy is this: Grandparent -> Parent -> Child
# The output would be:
# Child: Parent
# Parent: Grandparent
# Grandparent: .
def direct_parent(is_a):
    lst = is_a.split(';')
    ids = [i.split(' ! ')[0] for i in lst]
    return ';'.join(ids)
output['direct_parent_ids'] = output.is_a.apply(direct_parent)
print(output)

# Add "indirect" relationships
# A indirect relationship looks like this:
# If the hierarchy is this: Grandparent -> Parent -> Child
# The output would be:
# Child: Parent;Grandparent
# Parent: Grandparent
# Grandparent: .
def indirect_parent(parent_ids):
    if parent_ids=='.':
        return '.'
    parents = [i for i in parent_ids.split(';') if i!='.']
    
    # Get parents of parents
    for parent in parents:
        grandparents = [i for i in output[output.id==parent]['direct_parent_ids'].to_string(index = False, header = False).strip().split(';') if i!='.']
        new_grandparents = list(set([i for i in grandparents if i not in parents]))
        parents += new_grandparents
    
    parents.sort()
    return ';'.join(parents)
output['indirect_parent_ids'] = output.direct_parent_ids.apply(indirect_parent)
print(output)

output.to_csv('Analysis_files/5_term_relationships.csv', index = False)