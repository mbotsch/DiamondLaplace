# Diamond Laplace 🔷

This repository contains the implementation of the paper [Bunge, Botsch, Alexa, "The Diamond Laplace for Polygonal and Polyhedral Meshes", Symposium on Geometry Processing 2021](https://ls7-gv.cs.tu-dortmund.de/downloads/publications/2021/sgp21.pdf).

Since we use the [pmp-library](http://www.pmp-library.org/) and [OpenVolumeMesh](https://www.graphics.rwth-aachen.de/software/openvolumemesh/) as submodules, you have to clone the repository recursively:

    git clone --recursive git@github.com:mbotsch/polyLaplace.git

Configure and build:

    cd polyLaplace && mkdir build && cd build && cmake .. && make

This will automatically build our code and all dependencies. Finally, start the app with a polygon mesh or a polyhedral mesh:

    ./diamond_laplacian ../data/Surface/boar.obj
    ./diamond_laplacian ../data/Volume/bunny.mesh

Have fun!
