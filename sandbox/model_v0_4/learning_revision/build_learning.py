"""Build only E, preserving the A-D sources and PDFs."""
from pathlib import Path
import argparse
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'sandbox/model_v0_4'))
from build_documents import HEADER

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--tectonic');args=ap.parse_args()
    docs=ROOT/'docs/modelo_v0_4';build=docs/'latex';out=ROOT/'output/pdf/modelo_v0_4'
    md=docs/'E_cuaderno_de_aprendizaje_v0_4.md'
    source=md.read_text(encoding='utf-8')
    for path in re.findall(r'!\[[^\]]*\]\(([^)]+)\)',source):
        assert (docs/path).is_file(),path
    header=build/'header_learning.tex'
    header.write_text(HEADER.replace('Revisión 0.4','Modelo 0.4 | E-r02'),encoding='utf-8')
    pandoc=shutil.which('pandoc')
    if not pandoc:
        import pypandoc
        pandoc=pypandoc.get_pandoc_path()
        if not Path(pandoc).is_file():pandoc+='.exe'
    meta=re.match(r'\A---\n.*?\n---\n',source,re.S);assert meta
    source=source[:meta.end()]+'\n\\newpage\n\n'+source[meta.end():]
    tex=build/(md.stem+'.tex')
    cmd=[str(pandoc),'-f','markdown+tex_math_dollars+raw_tex','-t','latex','--standalone',
         '--toc','--toc-depth=1','--resource-path',str(docs),'-V','documentclass=article',
         '-V','fontsize=10pt','-V','papersize=a4','-V','geometry:top=22mm,bottom=22mm,left=22mm,right=22mm',
         '-V','colorlinks=true','-V','linestretch=1.07','-V','lang=es','--include-in-header',str(header),'-o',str(tex)]
    subprocess.run(cmd,input=source,text=True,encoding='utf-8',check=True,cwd=docs)
    content=tex.read_text(encoding='utf-8').replace(str(docs).replace('\\','/')+'/', '')
    content=content.replace('{figuras/', '{../figuras/')
    content=re.sub(r'\\texttt\{([^{}]*)\}',lambda m:r'\nolinkurl{'+re.sub(r'\\([_#$%&])',r'\1',m.group(1))+'}',content)
    tex.write_text(content,encoding='utf-8')
    if args.tectonic:
        subprocess.run([str(Path(args.tectonic).resolve()),str(tex),'--outdir',str(out),'--keep-logs'],check=True,cwd=build)
        shutil.copy2(out/(md.stem+'.log'),docs/'verificaciones/compilacion'/(md.stem+'.log'))
    print(tex)

if __name__=='__main__':main()
