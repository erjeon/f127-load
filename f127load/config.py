#!/usr/bin/env python3
"""Read and write the system file.

JSON was what the wizard wrote first. It is fine for a program and poor for a
person: it cannot carry a comment, so a number like a box edge sat there with no
way to say where it came from or what happens if it changes. Anyone editing the
file by hand had to go back to the wizard to find out.

The format here is one setting a line, with units left in and a comment column:

    box             24.03 nm     # 5.02 nm clear of the box face

Written by hand it is easier to change than JSON, and read back it is the same
dictionary. JSON is still accepted, so a config written by an earlier version
still builds.

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
