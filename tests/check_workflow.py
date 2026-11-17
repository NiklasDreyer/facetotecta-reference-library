#!/usr/bin/env python3
# =============================================================================
# check_workflow.py: static checks of the workflow (no tools needed).
#
#   1. every Python script compiles
#   2. every YAML file (config, profile, conda environments) parses
#   3. every script called from a rule exists, and every "--option" a rule
#      passes to a Python or R script is defined in that script
#   4. every {placeholder} in a rule's shell command refers to input, output,
#      params, threads, log, wildcards or a known variable, and every
#      {params.x}, {input.x} and {output.x} used is defined in that rule
#
# Run from the repository root:  python tests/check_workflow.py
# Exits with status 1 if a problem is found.
# =============================================================================

import glob
import os
import py_compile
import re
import sys

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WF = os.path.join(ROOT, "workflow")
problems = []

# 1. Python scripts compile.
for f in glob.glob(os.path.join(WF, "scripts", "*.py")) + glob.glob(os.path.join(ROOT, "tests", "*.py")):
    try:
        py_compile.compile(f, doraise=True)
    except py_compile.PyCompileError as e:
        problems.append(f"compile error: {e}")

# 2. YAML files parse.
for f in ([os.path.join(ROOT, "config", "config.yaml"), os.path.join(ROOT, "profiles", "slurm", "config.yaml")]
          + glob.glob(os.path.join(WF, "envs", "*.yaml"))):
    try:
        yaml.safe_load(open(f))
    except yaml.YAMLError as e:
        problems.append(f"YAML error in {f}: {e}")


def script_options(path):
    """Options defined by a Python (argparse) or R (optparse) script."""
    src = open(path).read()
    if path.endswith(".py"):
        return set(re.findall(r'add_argument\(\s*"(--[\w-]+)"', src))
    return set(re.findall(r'make_option\("(--[\w_]+)"', src))


SHELL = re.compile(r'shell:\s*\n((?:\s+(?:#.*|".*")\s*\n)+)')
ROOTS = {"input", "output", "params", "threads", "log", "wildcards", "SCRIPTS", "OUT", "resources", "PB"}
rule_files = glob.glob(os.path.join(WF, "rules", "*.smk"))

for smk in rule_files:
    txt = open(smk).read()
    name = os.path.relpath(smk, ROOT)
    for m in SHELL.finditer(txt):
        block = "".join(re.findall(r'"((?:[^"\\]|\\.)*)"', m.group(1)))
        # 3. Script calls and options.
        for cm in re.finditer(r'(?:python|Rscript|bash) \{SCRIPTS\}/([\w.]+)(.*?)(?=;|>|\|\||$)', block):
            script, args = cm.group(1), cm.group(2)
            path = os.path.join(WF, "scripts", script)
            if not os.path.exists(path):
                problems.append(f"{name}: missing script {script}")
                continue
            if script.endswith((".py", ".R")):
                for opt in set(re.findall(r'(--[\w-]+)', args)) - script_options(path):
                    problems.append(f"{name}: {script} does not define {opt}")
        # 4. Placeholders.
        flat = block.replace("{{", "").replace("}}", "")
        for ph in re.findall(r'\{([^{}]*)\}', flat):
            if re.split(r'[.\[:]', ph)[0] not in ROOTS:
                problems.append(f"{name}: unknown placeholder {{{ph}}}")

    # 4b. Names used in shell commands are defined in the same rule.
    for rm in re.finditer(r'\n(?:rule (\w+)|    rule):(.*?)(?=\nrule |\n    rule:|\Z)', txt, re.S):
        body = rm.group(2)
        for section in ("params", "output", "input"):
            sm = re.search(rf'\n\s+{section}:\n(.*?)(?=\n\s+(?:input|output|params|threads|log|conda|shell|wildcard_constraints|resources|default_target):|\Z)', body, re.S)
            defined = set(re.findall(r'^\s+(\w+)=', sm.group(1), re.M)) if sm else set()
            used = set(re.findall(rf'\{{{section}\.(\w+)', body))
            if section == "input" and "unpack(" in body:
                continue
            for u in used - defined:
                problems.append(f"{name}: rule {rm.group(1) or '(loop)'} uses undefined {section}.{u}")

if problems:
    print("\n".join(problems))
    sys.exit(1)
print(f"All checks passed ({len(rule_files)} rule files).")
