mkdir build
cd build
cmake ..
make
./manifold ../examples/input.obj ../examples/manifold.obj
#./simplify -i ../examples/manifold.obj -o ../examples/simplified_manifold.obj -c 1e-2 -f 10000 -r 0.2
cd ..

