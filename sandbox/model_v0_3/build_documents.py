"""Convert the five editable Markdown documents to LaTeX for Tectonic.

Usage: python build_documents.py [--compile]
Requires Pandoc (or pypandoc_binary) and Tectonic for --compile.
Run at any working directory; all paths are relative to the repository.
"""
from pathlib import Path
import argparse
import json
import re
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[2]
DOCS=ROOT/'docs/modelo_v0_3'
BUILD=DOCS/'latex'
PDF=ROOT/'output/pdf/modelo_v0_3'

HEADER=r'''
\usepackage{amsmath,amssymb,mathtools}
\usepackage{microtype}
\usepackage{fancyhdr}
\usepackage{needspace}
\usepackage{etoolbox}
\usepackage{float}
\usepackage{caption}
\captionsetup{font=small,labelformat=empty,skip=7pt}
\floatplacement{figure}{H}
\definecolor{DocBlue}{HTML}{173F57}
\definecolor{DocGray}{HTML}{52616C}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\sffamily\color{DocGray} pySNSPD | Desarrollo del modelo}
\fancyhead[R]{\small\sffamily\color{DocGray} Revisión 0.3}
\fancyfoot[L]{\small\sffamily\color{DocGray} 8 de septiembre de 2026}
\fancyfoot[R]{\small\sffamily\thepage}
\setlength{\headheight}{15pt}
\setlength{\emergencystretch}{2em}
\setlength{\parskip}{5pt plus 1pt}
\setlength{\parindent}{0pt}
\widowpenalty=10000
\clubpenalty=10000
\displaywidowpenalty=10000
\predisplaypenalty=10000
\renewcommand{\arraystretch}{1.16}
\AtBeginEnvironment{longtable}{\small}
\pretocmd{\section}{\Needspace{6\baselineskip}}{}{}
\pretocmd{\subsection}{\Needspace{5\baselineskip}}{}{}
\AtBeginDocument{\hypersetup{colorlinks=true,linkcolor=DocBlue,urlcolor=DocBlue}}
\allowdisplaybreaks[1]
'''

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--compile',action='store_true');args=ap.parse_args()
    BUILD.mkdir(parents=True,exist_ok=True);PDF.mkdir(parents=True,exist_ok=True)
    pandoc=shutil.which('pandoc')
    if not pandoc:
        import pypandoc
        pandoc=pypandoc.get_pandoc_path()
        if not Path(pandoc).exists() and Path(pandoc+'.exe').exists():pandoc+='.exe'
    header=BUILD/'header.tex';header.write_text(HEADER,encoding='utf-8')
    results=[]
    for md in sorted(DOCS.glob('[ABCDE]_*_v0_3.md')):
        source=md.read_text(encoding='utf-8')
        images=re.findall(r'!\[[^\]]*\]\(([^)]+)\)',source)
        for image in images:
            if not (DOCS/image).is_file():raise FileNotFoundError(image)
        # Preserve user-authored mathematics, require figures as physical files.
        meta=re.match(r'\A---\n.*?\n---\n',source,re.S)
        assert meta,md
        source=source[:meta.end()]+'\n\\newpage\n\n'+source[meta.end():]
        tex=BUILD/(md.stem+'.tex')
        cmd=[str(pandoc),'-f','markdown+tex_math_dollars+raw_tex','-t','latex','--standalone',
             '--toc','--toc-depth=2','--top-level-division=section','--resource-path',str(DOCS),
             '-V','documentclass=article','-V','fontsize=10pt','-V','papersize=a4',
             '-V','geometry:top=22mm,bottom=22mm,left=22mm,right=22mm',
             '-V','colorlinks=true','-V','linestretch=1.07','-V','lang=es',
             '--include-in-header',str(header),'-o',str(tex)]
        subprocess.run(cmd,input=source,text=True,encoding='utf-8',check=True,cwd=DOCS)
        # Pandoc resource paths can be absolute. LaTeX is shipped alongside figures.
        content=tex.read_text(encoding='utf-8').replace(str(DOCS).replace('\\','/')+'/', '')
        content=content.replace('{figuras/', '{../figuras/')
        if md.name.startswith('C_'):
            content=content.replace(r'\setcounter{tocdepth}{2}',r'\setcounter{tocdepth}{1}')
        # Long code paths and commit hashes must remain readable at line breaks.
        def breakable_code(match):
            value=re.sub(r'\\([_#$%&])',r'\1',match.group(1))
            return r'\nolinkurl{'+value+'}'
        content=re.sub(r'\\texttt\{([^{}]*)\}',breakable_code,content)
        tex.write_text(content,encoding='utf-8')
        if args.compile:
            subprocess.run(['tectonic',str(tex),'--outdir',str(PDF),'--keep-logs'],cwd=DOCS,check=True)
        results.append({'document':md.name,'latex':tex.name,'figures':images})
    (BUILD/'sources.json').write_text(json.dumps(results,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(results,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
