import re

def normalize_sample_id(name):
    name = name.replace(' (', '(')
    name = name.replace(' ', '-')
    name = re.sub('(?<!ELA1912|ELA1923|ELA1924)-?ADN[1-9]$', '', name)
    name = re.sub('-?ambBED$', '', name)
    name = re.sub('EMADN', 'EM', name)
    name = re.sub('(ELA|EM)-', r'\1', name)
    return name