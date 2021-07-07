# Spheres PDE solver

Preliminary version, no MPI support in the output.

## 

1. Checkout repo
```
git clone https://github.com/vovannikov/spheres.git
```

2. Create directory for build and get in there
```
mkdir spheres-build
cd spheres-build
```

3. Configure the application providing `DEAL_II_DIR` and `CMAKE_BUILD_TYPE` (Release or Debug)
```
cmake ../spheres -DDEAL_II_DIR=/path/to/dealii -DCMAKE_BUILD_TYPE=Release
```

4. Build the target
```
make
```

5. Run with given settings
```
./src/spheres_runner ../spheres/settings/settings.prm [2d or 3d]
```