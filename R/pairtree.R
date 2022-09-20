

# install --------
# ** deps -----
#' conda update -n base conda
#' conda create -n pairtree
#' conda install numpy
#' conda install scikit-learn scipy
#' conda install tqdm colorlover numba
#' conda install plotly # flor plotting
#' 
# ** clone and install -------
#' Clone the Pairtree repository, then download and build the 
#' C code required to fit subclone frequencies to the tree. 
#' This algorithm was published in Jia et al., and uses the 
#' authors' implementation with minor modifications.
#' 
#' cd ~/apps/pairtree
#' module load git
#' git clone https://github.com/jwintersinger/pairtree
#' cd pairtree/lib
#' git clone https://github.com/jwintersinger/projectppm
#' cd projectppm
#' bash make.sh
#' 
# ** test installation --------
#' PTDIR=$HOME/apps/pairtree/pairtree
#' cd $PTDIR/example
#' mkdir results && cd results
#' # Run Pairtree.
#' $PTDIR/bin/pairtree --params $PTDIR/example/example.params.json $PTDIR/example/example.ssm example.results.npz
# Plot results in an HTML file.
#' $PTDIR/bin/plottree --runid example $PTDIR/example/example.ssm $PTDIR/example/example.params.json example.results.npz example.results.html
#' # View the HTML file.
#' firefox example.results.html