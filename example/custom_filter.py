def custom_filter(record):
    info = record.INFO
    support = info['SU'][0]
    paired_end_support = info['PE'][0]
    split_read_support = info['SR'][0]
    if support > 10 and paired_end_support > 5 and split_read_support > 5:
        yield record 
