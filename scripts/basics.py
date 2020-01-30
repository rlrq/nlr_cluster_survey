def head(x, n=5):
    return tuple(x)[:n]

def tail(x, n=5):
    return tuple(x)[-n:]

## sorting key ('key' variable) should work on 'x' in 'for x in iterable'
def ties(iterable, sorting_function, tie_breaker = lambda x: x.items()[0], **kwargs):
    wanted_x = sorting_function(iterable, **kwargs)
    max_items = {k: v for k, v in iterable if key(k) == value} if type(iterable, dict) else \
                type(iterable)([x for x in iterables if sorting_key(x) == sorting_key(wanted_x)])
    return tie_breaker(max_items)

def print_iter(iterable):
    for entry in iterable:
        print(entry)
    return
