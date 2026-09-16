#!/usr/bin/env python3
"""Terminal dressing, kept in one place.

Colour is switched off when the output is not a terminal, so a log file or a
pipe stays readable. NO_COLOR is honoured. The palette matches the figures in
the accompanying paper: navy for the core, grey for the corona, blue for water,
red for the solute.
"""
import os
import shutil
import sys
import unicodedata

_ON = sys.stdout.isatty() and not os.environ.get("NO_COLOR")


def _c(code):
    return (lambda s: f"\033[{code}m{s}\033[0m") if _ON else (lambda s: s)


bold = _c("1")
dim = _c("38;5;245")          # corona grey
water = _c("38;5;68")         # water blue, the accent
solute = _c("38;5;167")       # solute red, for anything that needs a second look
good = _c("38;5;71")
warn = _c("38;5;179")
core = _c("38;5;60")          # core navy

WIDTH = min(shutil.get_terminal_size((80, 24)).columns, 76)


def w(s):
    """Display width. Hangul and CJK take two columns, so len() misaligns boxes."""
    return sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in str(s))


def pad(s, n):
    return str(s) + " " * max(0, n - w(s))


def rule(char="─"):
    print(dim(char * WIDTH))


def banner(title, subtitle=""):
    print()
    print(core("  ┌" + "─" * (WIDTH - 4) + "┐"))
    print(core("  │ ") + bold(water(pad(title, WIDTH - 6))) + core(" │"))
    if subtitle:
        print(core("  │ ") + dim(pad(subtitle, WIDTH - 6)) + core(" │"))
    print(core("  └" + "─" * (WIDTH - 4) + "┘"))


def step(n, total, title):
    print()
    print(f"  {water(f'[{n}/{total}]')} {bold(title)}")
    print(f"  {dim('─' * (WIDTH - 4))}")


def note(text, kind="info"):
    mark = {"info": dim("·"), "ok": good("✓"), "warn": warn("!"),
            "bad": solute("×")}[kind]
    body = {"info": dim, "ok": good, "warn": warn, "bad": solute}[kind]
    print(f"    {mark} {body(text)}")


def value(label, val, unit="", tone=None):
    """One aligned label and value, for a result the user should register."""
    paint = {"good": good, "warn": warn, "bad": solute, None: bold}[tone]
    print(f"    {dim(pad(label, 26))} {paint(str(val))}{dim(' ' + unit if unit else '')}")


def table(rows, headers=None):
    cols = len(rows[0]) if rows else 0
    wid = [0] * cols
    body = ([headers] if headers else []) + [[str(c) for c in r] for r in rows]
    for r in body:
        for i, c in enumerate(r):
            wid[i] = max(wid[i], w(c))
    if headers:
        print("    " + dim("  ".join(pad(h, wid[i]) for i, h in enumerate(headers))))
        print("    " + dim("  ".join("─" * x for x in wid)))
    for r in rows:
        print("    " + "  ".join(pad(c, wid[i]) for i, c in enumerate(r)))


def prompt(text, default=None):
    tail = dim(f" [{default}]") if default is not None else ""
    v = input(f"    {water('›')} {text}{tail} ").strip()
    if not sys.stdin.isatty():
        print()          # a pipe echoes no newline, which runs the log together
    return v


def choices(options):
    for i, (_, label) in enumerate(options, 1):
        print(f"      {water(str(i) + ')')} {label}")
