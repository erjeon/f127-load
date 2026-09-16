#!/usr/bin/env python3
"""Find the config the user most likely means, and say which one it is.

The browser page hands you a file through the download folder. Then you have to
work out where it landed, change to the directory the tool is in, and type a path
that joins the two. That is three steps between deciding what to build and
building it, and none of them is about the simulation.

With no argument, look in the current directory first, then in the usual
download folders, and take the newest file that actually parses as one of these
configs. Reading the file is what makes this safe: a JSON that is not an f127
config is skipped rather than half-built.

    python scripts/find_config.py            the path, on stdout
    python scripts/find_config.py --explain  the path and where it came from

Pukyong National University / NCHM Lab.  Eunryul Jeon <qlsguswjs@pukyong.ac.kr>
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from f127load import config

# in order of how likely it is that this is the one meant
PLACES = [Path.cwd(), Path.home() / "Downloads", Path.home() / "Desktop",
          Path.home() / "다운로드"]
NAMES = ("system.cfg", "system.json")
REQUIRED = ("route",)


def looks_like_one(p):
    try:
        c = config.read(p)
    except Exception:
        return False
    if not isinstance(c, dict):
        return False
    # box and chains under either spelling, and a solute, and where it starts
    has_box = any(k in c for k in ("box", "box_nm"))
    has_solute = any(k in c for k in ("solute", "guest"))
    return has_box and has_solute and all(k in c for k in REQUIRED)


def ago(p):
    s = time.time() - p.stat().st_mtime
    if s < 90:
        return "just now"
    if s < 5400:
        return f"{round(s / 60)} minutes ago"
    if s < 172800:
        return f"{round(s / 3600)} hours ago"
    return f"{round(s / 86400)} days ago"


def find():
    # the conventional names where you are standing win, whatever their age
    for name in NAMES:
        p = Path.cwd() / name
        if p.is_file() and looks_like_one(p):
            return p, "in this directory"
    best, where = None, ""
    for place in PLACES:
        if not place.is_dir():
            continue
        for p in place.glob("*.json"):
            if p.is_file() and looks_like_one(p) and (best is None
                                                      or p.stat().st_mtime > best.stat().st_mtime):
                best, where = p, place.name
        for p in place.glob("*.cfg"):
            if p.is_file() and looks_like_one(p) and (best is None
                                                      or p.stat().st_mtime > best.stat().st_mtime):
                best, where = p, place.name
    if best is None:
        return None, ""
    return best, f"in {where}, saved {ago(best)}"


def main():
    p, where = find()
    if p is None:
        return 1
    if "--explain" in sys.argv:
        print(f"{p}\t{where}")
    else:
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
