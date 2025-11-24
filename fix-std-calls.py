#!/usr/bin/env python3
from __future__ import annotations
import sys, re, argparse
from pathlib import Path
from dataclasses import dataclass
from collections import defaultdict

DEFAULT_FUNCS = {"min", "max", "log", "pow", "cbrt", "atan", "sinh", "cosh", "abs", "fabs", "sin", "cos", "tan", "asin", "acos", "sqrt", "exp", "tanh"}

ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")

ERROR_WITH_FILE_RE = re.compile(
    r"""^
        (?P<file>.+?\.(?:C|cc|cpp|cxx))   # C++ file
        :(?P<line>\d+)                    # line
        (?::\d+)?                         # optional column
        :\s*error:\s*
        .*?no\ matching\ function\ for\ call\ to\s*'(?P<func>\w+)\s*\(
    """,
    re.VERBOSE,
)
CARRIER_RE = re.compile(r"^(?P<file>.+?\.(?:C|cc|cpp|cxx)):(?P<line>\d+)(?::\d+)?:")
NO_MATCH_FUNC_RE = re.compile(r"no\ matching\ function\ for\ call\ to\s*'(?P<func>\w+)\s*\(", re.IGNORECASE)

KW_BLOCK_HEADERS = {"if", "for", "while", "switch", "catch"}
FUNC_SIG_RE = re.compile(
    r"""
    ^(?:[ \t]*(?:template\s*<[^>]*>\s*)*)?
    (?:[\w:\<\>\[\]\s\*&]+?)?
    \b([A-Za-z_]\w*(?:::[A-Za-z_]\w*)*)
    \s*\([^;{}]*\)
    (?:\s*(?:const|noexcept|constexpr|override|final|->\s*[\w:\s<>\*&]+))*\s*
    \{
    """,
    re.VERBOSE | re.MULTILINE,
)

@dataclass
class ErrorHit:
    file: Path
    line: int
    func: str

def strip_line_comments(s: str) -> str:
    i = s.find("//")
    return s if i < 0 else s[:i]

def find_enclosing_block_bounds(src: str, err_line_idx0: int) -> tuple[int, int]:
    lines = src.splitlines(True)
    upto = sum(len(x) for x in lines[:err_line_idx0 + 1])
    depth = 0
    opens = []
    pos = 0
    for raw in lines:
        line = strip_line_comments(raw)
        for j, ch in enumerate(line):
            if ch == "{":
                opens.append(pos + j); depth += 1
            elif ch == "}":
                if depth:
                    opens.pop(); depth -= 1
        pos += len(raw)
        if pos >= upto:
            break
    if not opens:
        return (0, len(src))
    open_pos = opens[-1]
    depth2, i, n = 1, open_pos + 1, len(src)
    while i < n and depth2 > 0:
        c = src[i]
        if c == "{": depth2 += 1
        elif c == "}":
            depth2 -= 1
            if depth2 == 0: return (open_pos, i)
        i += 1
    return (0, len(src))

def looks_like_function_header(header: str) -> bool:
    t = header.strip()
    if ";" in t or "(" not in t or ")" not in t: return False
    head = t.split("(")[0].split()[-1] if "(" in t else ""
    return head not in KW_BLOCK_HEADERS

def find_function_block(src: str, err_line_idx0: int) -> tuple[int, int]:
    lines = src.splitlines(True)
    open_idx, close_idx = find_enclosing_block_bounds(src, err_line_idx0)
    acc = 0; char_line = []
    for i, ln in enumerate(lines):
        char_line.append((acc, i)); acc += len(ln)
    open_line = 0
    for cp, ln in reversed(char_line):
        if cp <= open_idx: open_line = ln; break
    start_ctx = max(0, open_line - 8)
    header_text = src[sum(len(x) for x in lines[:start_ctx]) : open_idx + 1]
    if FUNC_SIG_RE.search(header_text) or looks_like_function_header(header_text):
        return open_idx, close_idx
    # one level outward
    outer_line = max(0, open_line - 1)
    open2, close2 = find_enclosing_block_bounds(src, outer_line)
    if (open2, close2) != (open_idx, close_idx):
        start_ctx2 = max(0, outer_line - 8)
        header2 = src[sum(len(x) for x in lines[:start_ctx2]) : open2 + 1]
        if FUNC_SIG_RE.search(header2) or looks_like_function_header(header2):
            return open2, close2
    return open_idx, close_idx

def already_has_using(block_src: str, name: str) -> bool:
    return re.search(rf"\busing\s+std::{re.escape(name)}\s*;", block_src) is not None

def add_using_decls(block_src: str, needed: set[str]) -> str:
    if not needed: return block_src
    after = block_src[1:]
    m = re.search(r"\n([ \t]*)\S", after)
    indent = m.group(1) if m else ""
    inject = "".join(f"\n{indent}using std::{n};" for n in sorted(needed))
    return block_src[:1] + inject + block_src[1:]

def rewrite_std_calls(block_src: str, names: set[str]) -> str:
    for n in sorted(names, key=len, reverse=True):
        block_src = re.sub(rf"\bstd::{re.escape(n)}\s*\(", f"{n}(", block_src)
    return block_src

def read_text_from(path_or_dash: str) -> str:
    if path_or_dash == "-":
        return sys.stdin.read()
    return Path(path_or_dash).read_text(errors="ignore")

def parse_log(text: str, targets: set[str], verbose: bool, trace: bool) -> tuple[list[ErrorHit], dict]:
    stats = dict(lines=0, carriers=0, direct_matches=0, fallback_matches=0)
    text = ANSI_RE.sub("", text)
    lines = text.splitlines()
    stats["lines"] = len(lines)
    hits: list[ErrorHit] = []

    # direct matches
    for ln in lines:
        m = ERROR_WITH_FILE_RE.match(ln)
        if m:
            func = m.group("func")
            if func in targets:
                fpath = Path(m.group("file")); line_no = int(m.group("line"))
                hits.append(ErrorHit(fpath, line_no, func))
                stats["direct_matches"] += 1
                if verbose: print(f"[match:direct] {fpath}:{line_no} func={func}", flush=True)
            elif trace:
                print(f"[trace:direct-skip] func={m.group('func')} not in targets", flush=True)

    if hits:
        return hits, stats

    # fallback association: carrier → next "no matching function" line
    last_file, last_line = None, None
    for ln in lines:
        cm = CARRIER_RE.match(ln)
        if cm:
            stats["carriers"] += 1
            last_file = Path(cm.group("file"))
            last_line = int(cm.group("line"))
            if trace: print(f"[trace:carrier] {last_file}:{last_line}", flush=True)
            continue
        nm = NO_MATCH_FUNC_RE.search(ln)
        if nm and last_file and last_line is not None:
            func = nm.group("func")
            if func in targets:
                hits.append(ErrorHit(last_file, last_line, func))
                stats["fallback_matches"] += 1
                if verbose: print(f"[match:fallback] {last_file}:{last_line} func={func}", flush=True)
            elif trace:
                print(f"[trace:fallback-skip] func={func} not in targets", flush=True)

    return hits, stats

def patch_file(path: Path, file_hits: list[ErrorHit], targets: set[str], verbose: bool, make_backup: bool) -> bool:
    if not path.exists():
        print(f"[warn] Missing file: {path}", flush=True)
        return False
    original = path.read_text(errors="ignore")
    src = original
    modified = False

    for h in sorted(file_hits, key=lambda x: x.line, reverse=True):
        err_line0 = max(0, h.line - 1)
        open_idx, close_idx = find_function_block(src, err_line0)
        block = src[open_idx:close_idx + 1]
        present_std_calls = {m.group(1) for m in re.finditer(r"\bstd::(\w+)\s*\(", block)}
        needed = ({h.func} | (present_std_calls & targets))
        to_add = {n for n in needed if not already_has_using(block, n)}
        new_block = add_using_decls(block, to_add)
        new_block = rewrite_std_calls(new_block, needed)
        if new_block != block:
            src = src[:open_idx] + new_block + src[close_idx + 1:]
            modified = True
            if verbose:
                print(f"[patched] {path} ~line {h.line}: using {sorted(to_add)}; rewrote {sorted(needed)}", flush=True)
        elif verbose:
            print(f"[noop] {path} ~line {h.line}: nothing to change", flush=True)

    if modified:
        if make_backup:
            bak = path.with_suffix(path.suffix + ".bak")
            if not bak.exists():
                bak.write_text(original)
        path.write_text(src)
        print(f"[ok] Patched {path}", flush=True)
    else:
        print(f"[skip] No changes needed for {path}", flush=True)
    return modified

def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("log", help="compiler log file path or '-' for stdin")
    ap.add_argument("--funcs", default=",".join(sorted(DEFAULT_FUNCS)),
                    help="comma-separated function names to fix")
    ap.add_argument("--verbose", action="store_true", help="print matches and patch actions")
    ap.add_argument("--trace", action="store_true", help="very chatty: show why things were skipped")
    ap.add_argument("--backup", action="store_true", help="write .bak files (off by default)")
    args = ap.parse_args(argv[1:])

    try:
        text = read_text_from(args.log)
    except Exception as e:
        print(f"[error] failed reading log '{args.log}': {e}", flush=True)
        return 2

    targets = {s.strip() for s in args.funcs.split(",") if s.strip()}
    hits, stats = parse_log(text, targets, args.verbose, args.trace)

    # group by file
    by_file = defaultdict(list)
    for h in hits:
        by_file[h.file].append(h)

    any_mod = False
    for fpath, flist in by_file.items():
        try:
            any_mod |= patch_file(fpath, flist, targets, args.verbose, make_backup=args.backup)
        except Exception as e:
            print(f"[error] Failed patching {fpath}: {e}", flush=True)

    # Always print a summary
    print(
        f"[summary] lines={stats['lines']} carriers={stats['carriers']} "
        f"direct_matches={stats['direct_matches']} fallback_matches={stats['fallback_matches']} "
        f"files={len(by_file)} patched={'yes' if any_mod else 'no'}",
        flush=True,
    )
    return 0

if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
