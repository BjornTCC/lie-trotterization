from copy import deepcopy

def add_dicts(dct1: dict, dct2: dict, inplace: bool = False) -> dict | None:
    if inplace:
        _res = dct1
    else:
        _res = deepcopy(dct1)

    for key, val in dct2.items():
        if key in _res.keys():
            _res[key] += val
        else:
            _res[key] = val

    if not inplace:
        return _res