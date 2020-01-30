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
    
def keep_best_hit_per_query_subject_pair(fname, fout, bs_col, q_col = 0, s_col = 1, delim = '\t', header = True, length_col = 3, min_length = 0):
    with open(fname, 'r') as f:
        data_raw = [x[:-1].split(delim) for x in f.readlines()]
    output = [delim.join(data_raw[0])] if header else []
    data = data_raw[1:] if header else data_raw
    data_rows = len(data)
    data.sort()
    last_row = data[0]
    last_q, last_s, last_bs = last_row[q_col], last_row[s_col], float(last_row[bs_col])
    for i,row in enumerate(data[1:]):
        if int(row[length_col]) < min_length:
            continue
        curr_q, curr_s, curr_bs = row[q_col], row[s_col], float(row[bs_col])
        if curr_q != last_q or curr_s != last_s:
            output.append(delim.join(last_row))
            last_row = row
            last_q, last_s, last_bs = last_row[q_col], last_row[s_col], float(last_row[bs_col])
        else:
            if curr_bs > last_bs:
                last_row = row
                last_q, last_s, last_bs = last_row[q_col], last_row[s_col], float(last_row[bs_col])
        if i == data_rows - 2:
            output.append(delim.join(last_row))
    f = open(fout, "w+")
    f.write('\n'.join(output))
    f.close()

def keep_query_aligned_to_selected_subjects(fname, fout, subjects, q_col = 0, s_col = 1, delim = '\t', header = True):
    with open(fname, 'r') as f:
        data_raw = [x[:-1].split(delim) for x in f.readlines()]
    output = [delim.join(data_raw[0])] if header else []
    data = data_raw[1:] if header else data_raw
    selected_queries = []
    for subject in subjects:
        selected_queries.extend([x[q_col] for x in data if subject in x[s_col]])
    output.extend([delim.join(x) for x in data if x[q_col] in selected_queries])
    f = open(fout, "w+")
    f.write('\n'.join(output))
    f.close()
