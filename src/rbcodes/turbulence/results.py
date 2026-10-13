"""
Results bookkeeping: numbers for a paper flow through results.json -> results.tex.

    from rbcodes.turbulence.results import Results, pm, asym, num, savefig
    R = Results('paper/results.json', 'paper/results.tex')
    R.update({'alpha-plane': pm(1.08, 0.46), 'rgeo': num(7.1)}, section='nb03')

In LaTeX: ``\\input{results.tex}`` and ``\\newcommand{\\R}[1]{\\csname R@#1\\endcsname}``; then
``\\R{alpha-plane}`` prints the value. Each update REPLACES its section, so a rerun never leaves
stale numbers. Ported from vsf_pipeline (``update_results``, ``pm``, ``asym``, ``num``, ``savefig``).
"""
import json
import os

__all__ = ['Results', 'pm', 'asym', 'num', 'savefig']


class Results:
    """results.json (sections of {key: LaTeX string}) and the generated results.tex."""

    def __init__(self, json_path='paper/results.json', tex_path='paper/results.tex'):
        self.json_path, self.tex_path = os.fspath(json_path), os.fspath(tex_path)

    def load(self):
        if not os.path.exists(self.json_path):
            return {}
        with open(self.json_path) as fh:
            allv = json.load(fh)
        if allv and not all(isinstance(v, dict) for v in allv.values()):
            allv = {'misc': {k: v for k, v in allv.items() if not isinstance(v, dict)}}
        return allv

    def update(self, values, section='misc'):
        """Store {key: value} under ``section`` (replacing that section) and rewrite results.tex."""
        allv = self.load()
        allv[section] = {k: str(v) for k, v in values.items()}
        os.makedirs(os.path.dirname(self.json_path) or '.', exist_ok=True)
        with open(self.json_path, 'w') as fh:
            json.dump(allv, fh, indent=1, sort_keys=True)
        self.write_tex(allv)
        return allv

    def flat(self, allv=None):
        allv = self.load() if allv is None else allv
        out = {}
        for sec in sorted(allv):
            out.update(allv[sec])
        return out

    def write_tex(self, allv=None):
        flat = self.flat(allv)
        os.makedirs(os.path.dirname(self.tex_path) or '.', exist_ok=True)
        with open(self.tex_path, 'w') as f:
            f.write('% auto-generated -- do not edit by hand\n')
            for k in sorted(flat):
                f.write('\\expandafter\\def\\csname R@%s\\endcsname{%s}\n' % (k, flat[k]))

    def __getitem__(self, key):
        return self.flat()[key]


def pm(x, e, nd=2):
    """'x ± e' in math mode."""
    return '\\ensuremath{%.*f \\pm %.*f}' % (nd, x, nd, e)


def asym(md, lo, hi, nd=2):
    """'md^{+(hi-md)}_{-(md-lo)}' in math mode."""
    return '\\ensuremath{%.*f^{+%.*f}_{-%.*f}}' % (nd, md, nd, hi - md, nd, md - lo)


def num(x, nd=1):
    return '\\ensuremath{%.*f}' % (nd, x)


def savefig(fig, name, outdir='paper/figures', dpi=150):
    """Save ``fig`` as outdir/name.pdf and .png."""
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, name + '.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(outdir, name + '.png'), bbox_inches='tight', dpi=dpi)
