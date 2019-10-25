mkdir AlphaPeel

# Assumes that the program and manual have both been built.

# To build the program run:
#cmake . ; make

# to build the manual:
#( cd docs ; make )

cp example_multi_locus_spec.txt AlphaPeel
cp example_single_locus_spec.txt AlphaPeel

cp -r example AlphaPeel

# Copy in the documentation.
cp docs/build/latex/AlphaPeel.pdf AlphaPeel

# Copy in the binaries
cp binaries/* AlphaPeel

# Create a version file

echo Version: 

zip -r AlphaPeel.zip AlphaPeel
