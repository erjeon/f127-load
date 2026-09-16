#!/usr/bin/env python3
"""Read and write the system file.

The system file is JSON, the same shape the browser designer downloads, so a
config from either source builds the same way. The wizard puts what it would
have said in a comment column into a "notes" object instead. The older one
setting per line format (system.cfg) is still read, so old configs still build.

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import json
import re
from pathlib import Path

UNITS = ("wt%", "nm", "mM", "M", "%", "ns")


def _value(raw):
    """Strip a trailing unit and give back a number where the text is one."""
    s = raw.strip()
    for u in sorted(UNITS, key=len, reverse=True):
        if s.endswith(" " + u) or s.endswith(u) and not s[:-len(u)].strip().isalpha():
            s = s[: len(s) - len(u)].strip()
            break
    if re.fullmatch(r"-?\d+", s):
        return int(s)
    if re.fullmatch(r"-?\d*\.\d+([eE][-+]?\d+)?", s):
        return float(s)
    return s


def read(path):
    """A settings file, or a JSON file written by an earlier version."""
    text = Path(path).read_text(encoding="utf-8")
    if text.lstrip().startswith("{"):
        return json.loads(text)
    out = {}
    for line in text.splitlines():
        line = line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        key, _, raw = line.partition(" ")
        if not raw.strip():
            continue
        out[key.strip()] = _value(raw)
    return out


def write_json(path, data):
    """The designer's shape, two space indent, one trailing newline."""
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                          encoding="utf-8")


def rows_to_json(rows, salts=()):
    """Turn the wizard's (key, value, unit, comment) rows into the designer's keys."""
    keymap = {"chains": "n_chains", "box": "box_nm", "guest": "solute",
              "count": "n_solute", "polymer": "wt_percent"}
    out, notes, salt = {}, {}, {}
    for key, value, unit, comment in rows:
        if not key:
            continue
        v = _value(f"{value} {unit}".strip())
        if key in salts:
            salt[key] = float(v)
        elif key == "salt" and value == "none":
            pass
        else:
            out[keymap.get(key, key)] = v
        if comment:
            notes[key] = comment
    out["salts"] = salt
    out["notes"] = notes
    return out


def write(path, rows, header=""):
    """rows: (key, value, unit, comment). A blank key leaves an empty line."""
    width = max((len(k) for k, *_ in rows if k), default=10)
    lines = []
    if header:
        lines += ["# " + h if h else "#" for h in header.splitlines()] + [""]
    for key, value, unit, comment in rows:
        if not key:
            lines.append("" if not comment else f"# {comment}")
            continue
        shown = f"{value} {unit}".strip()
        left = f"{key:<{width}}  {shown}"
        lines.append(f"{left:<34}# {comment}" if comment else left)
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
