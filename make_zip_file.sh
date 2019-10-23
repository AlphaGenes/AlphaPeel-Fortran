mkdir AlphaPeel

# Assumes that the program and manual have both been built.

# To build the program run:
#cmake . ; make

# to build the manual:
#( cd docs ; make )

cp example_multi_locus_spec.txt AlphaPeel
cp example_single_locus_spec.txt AlphaPeel

cp -r example AlphaPeel

# Need to figure out how to do the zip files.
cp docs/build/latex/AlphaPeel.pdf AlphaPeel

zip -r AlphaPeel.zip AlphaPeel
