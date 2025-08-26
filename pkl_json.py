#!/usr/bin/env python3
import pickle, json
from fractions import Fraction
import sys
import numpy as np
import sympy
from sympy import Rational as SympyRational, Integer as SympyInteger, Float as SympyFloat

PKL_PATH  = "./profile_dicts/profile_6.pkl"
JSON_PATH = "./profile_dicts/profile_6.json"

def to_jsonable(obj, _path="root"):
    
    """
    Recursively convert Python objects into JSON‐serializable primitives,
    handling built-in, numpy scalars, SymPy numbers, Fraction, and nested structures.
    Logs fallback cases to stderr for debugging.
    """

    # 1. Bool first
    if isinstance(obj, bool):
        return obj

    if isinstance(obj, np.ndarray):
        if obj.ndim == 0:
            scalar = obj.item()
            return to_jsonable(scalar, _path=_path + ".ndarray0d")
        else:
            # Example: convert array to nested lists; if you don't need arrays, you can skip or handle differently
            return [to_jsonable(e, _path=_path + ".ndarray_elem") for e in obj.tolist()]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)

    if isinstance(obj, SympyRational):
        # as_numer_denom returns (num, den) as SympyIntegers
        num, den = obj.as_numer_denom()
        return {"__fraction__": [int(num), int(den)]}
    # SymPy Integer (includes Sympy.One, Sympy.Zero, etc.)
    if isinstance(obj, SympyInteger):
        return int(obj)
    # SymPy Float
    if isinstance(obj, SympyFloat):
        return float(obj)

    if isinstance(obj, Fraction):
        return {"__fraction__": [obj.numerator, obj.denominator]}

    # 5. dict
    if isinstance(obj, dict):
        out = {"__dict__": []}
        for k, v in obj.items():
            k_json = to_jsonable(k, _path=_path + f".dict_key({repr(k)})")
            v_json = to_jsonable(v, _path=_path + f".dict_val({repr(k)})")
            out["__dict__"].append([k_json, v_json])
        return out

    # 6. tuple
    if isinstance(obj, tuple):
        return {"__tuple__": [to_jsonable(e, _path=_path + f".tuple[{i}]") for i, e in enumerate(obj)]}

    # 7. set
    if isinstance(obj, set):
        return {"__set__": [to_jsonable(e, _path=_path + f".set_elem({repr(e)})") for e in obj]}

    # 8. list
    if isinstance(obj, list):
        return [to_jsonable(e, _path=_path + f".list[{i}]") for i, e in enumerate(obj)]

    # 9. Built-in primitives
    if isinstance(obj, (int, float, str)) or obj is None:
        return obj

    # 10. Fallback: log for debugging
    print(f"Warning: fallback at {_path}: type={type(obj)}, repr={repr(obj)}", file=sys.stderr)
    return repr(obj)

def gen_json():
    # 1) Load pickle
    with open(PKL_PATH, "rb") as f:
        data = pickle.load(f)

    # 2) Convert recursively
    jsonable = to_jsonable(data)

    # 3) Dump JSON
    with open(JSON_PATH, "w") as f:
        json.dump(jsonable, f, indent=1)

    print(f"JSON written to {JSON_PATH}")

if __name__ == "__main__":
    gen_json()
