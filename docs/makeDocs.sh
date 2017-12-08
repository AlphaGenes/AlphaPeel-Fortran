docutils/docutils/tools/rst2latex.py --stylesheet test.sty AlphaPeelDocs.rst AlphaPeelDocs.tex 
pdflatex AlphaPeelDocs.tex
bibtex AlphaPeelDocs.tex
makeindexAlphaPeelDocs.tex
pdflatex AlphaPeelDocs.tex
pdflatex AlphaPeelDocs.tex

open AlphaPeelDocs.pdf