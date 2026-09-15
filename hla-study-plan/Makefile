.PHONY: all clean
all: plan.pdf
plan.pdf: plan.tex VERSION
	latexmk -pdf -interaction=nonstopmode -halt-on-error plan.tex
clean:
	latexmk -c plan.tex
