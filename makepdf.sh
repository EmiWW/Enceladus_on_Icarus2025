#!/bin/bash
pdflatex enceladus_crater.tex
bibtex enceladus_crater
pdflatex enceladus_crater.tex
pdflatex enceladus_crater.tex
