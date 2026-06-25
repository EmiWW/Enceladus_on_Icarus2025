#!/bin/bash

if [ "$1" = "full" ]; then
  pdflatex enceladus_crater.tex
  bibtex enceladus_crater
  pdflatex enceladus_crater.tex
  pdflatex enceladus_crater.tex
  cp enceladus_crater.pdf latest_copy.pdf
elif [ "$1" = "diff" ]; then
  latexdiff first_submission.tex enceladus_crater.tex >diff.tex
  pdflatex diff.tex
  bibtex diff
  pdflatex diff.tex
  pdflatex diff.tex
else
  pdflatex enceladus_crater.tex
  cp enceladus_crater.pdf latest_copy.pdf
fi
