#!/bin/bash

if [ "$1" = "full" ]; then
  pdflatex enceladus_crater.tex
  bibtex enceladus_crater
  pdflatex enceladus_crater.tex
  pdflatex enceladus_crater.tex
  cp enceladus_crater.pdf latest_copy.pdf
elif [ "$1" = "diff" ]; then
  latexdiff first_submission.tex enceladus_crater.tex >diff.tex
  #latexdiff --append-safecmd=url first_submission_edit.tex enceladus_crater.tex >diff.tex
  #latexdiff --append-textcmd="DIFadd" first_submission_edit.tex enceladus_crater.tex >diff.tex
  echo "diff.text has been made"
  pdflatex diff.tex
  bibtex diff
  pdflatex diff.tex
  pdflatex diff.tex
  cp diff.pdf latest_diff.pdf
else
  pdflatex enceladus_crater.tex
  cp enceladus_crater.pdf latest_copy.pdf
fi
