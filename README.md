# Watertight Manifold

Source code for the paper:

Huang, Jingwei, Hao Su, and Leonidas Guibas. [**Robust Watertight Manifold Surface Generation Method for ShapeNet Models.**](https://arxiv.org/abs/1802.01698) , arXiv preprint arXiv:1802.01698 (2018).

## Install and Run

For linux and mac user, run **sh demo.sh** to build and try the manifold example.

### Install

git clone --recursive -j8 git://github.com/foo/bar.git

mkdir build

cd build

cmake ..

make

### Manifold Software

We take a triangle mesh "input.obj" and generate a manifold "output.obj". The resolution is the number of leaf nodes of octree. The face number increases linearly with the resolution.

./manifold input.obj output.obj [resolution (Default 20000)]

### Simplify Algorithm

Our manifold software generates uniform manifold. For efficiency purpose, a mesh simplification can be used.

./simplify -i input.obj -o output.obj [-m] [-f face_num] [-c max_cost] [-r max_ratio]

-m 					// turn on manifold check, we don't output model if check fails

-f face_num 		// add termination condition when current_face_num <= face_num

-c max_cost			// add termination condition when quadric error >= max_cost

-r max_ratio		// add termination condition when current_face_num / origin_face_num <= max_ratio

### Example:

./simplify -i input.obj -o output.obj -m -c 1e-2 -f 10000 -r 0.2


## Authors
- [Jingwei Huang](mailto:jingweih@stanford.edu)

&copy; Jingwei Huang, Stanford University


**IMPORTANT**: If you use this software please cite the following in any resulting publication:
```
@article{huang2018robust,
  title={Robust Watertight Manifold Surface Generation Method for ShapeNet Models},
  author={Huang, Jingwei and Su, Hao and Guibas, Leonidas},
  journal={arXiv preprint arXiv:1802.01698},
  year={2018}
}
```
