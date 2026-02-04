#!/bin/bash

if ["$1" = "full" ]; then
  pdflatex enceladus_crater.tex
  bibtex enceladus_crater
  pdflatex enceladus_crater.tex
  pdflatex enceladus_crater.tex
else
  pdflatex enceladus_crater.tex
fi

#pdflatex new.tex
#bibtex enceladus_crater
#pdflatex new.tex
#pdflatex new.tex
