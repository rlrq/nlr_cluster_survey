import sys
import fileinput

def filter_column(fname, fitems, col, fout, inverse = False, delim = '\t'):
    with open(fname, 'r') as f:
        dat = [x[:-1].split(delim) for x in f.readlines()]
    items = []
    with open(fitems, 'r') as f:
        items_raw = [x[:-1].split(delim) for x in f.readlines()]
        for item in items_raw:
            items.extend(item)
    output = [delim.join(x) for x in dat if x[col] in items] if not inverse else \
             [delim.join(x) for x in dat if x[col] not in items]
    f = open(fout, "w+")
    f.write('\n'.join(output))
    f.close()

def single_col_to_list(fname, delim = '\t'):
    """
    reads a single-column file and generates a list of each line (removes trailing newline char)
    """
    with open(fname, 'r') as f:
        data = [x[:-1] for x in f.readlines()]
    return data

def cat_files(fnames, fout):
    with open(fout, "w+") as f, fileinput.input(fnames) as fin:
        for line in fin:
            f.write(line)
    return fout
    
