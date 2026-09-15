.PHONY: all clean
all: plan.pdf technical-notes.pdf
plan.pdf: plan.tex VERSION
	latexmk -pdf -interaction=nonstopmode -halt-on-error plan.tex
technical-notes.pdf: technical-notes.tex VERSION
	latexmk -pdf -interaction=nonstopmode -halt-on-error technical-notes.tex
clean:
	latexmk -c plan.tex technical-notes.tex
