# Bone PDE solver

Preliminary version, no MPI support in the output.

## 

1. Checkout repo
```
git clone https://github.com/vovannikov/bone.git
```

2. Create directory for build and get in there
```
mkdir bone-build
cd bone-build
```

3. Configure the application providing `DEAL_II_DIR` and `CMAKE_BUILD_TYPE` (Release or Debug)
```
cmake ../bone -DDEAL_II_DIR=/path/to/dealii -DCMAKE_BUILD_TYPE=Release
```

4. Run the example
```
./src/bone_example
```

If you want to save the VTK output for debug purpose, then provide the folder for results:
```
./src/bone_example /path/to/save/results
```