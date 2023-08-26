## One-Dimensional Helmholtz Equation in the Cartesian Coordinate System

We solve the following equation using the Finite Difference Method(FDM)

$$\frac{d^2 u}{dx^2} - \lambda u = f,\quad \lambda \geq 0$$

with boundary conditions 
$$\alpha_l \frac{du}{dn}(x_l) + \beta_l u(x_l)=g_l$$
and
$$\alpha_r \frac{du}{dn}(x_r) + \beta_r u(x_r)=g_r$$

Note that $\frac{du}{dn}(x_l) = -\frac{du}{dx}(x_l)$, and $\frac{du}{dn}(x_r) = \frac{du}{dx}(x_r)$.

To build and run the code, invoke the following commands
```
make TARGET=TestNonPeriodic1D_MPI.cpp
make run Nproc=4
```
Here, Nproc denotes the number of processors used to run the code. To use the periodic boundary condition along the x-direction, build and run the below code

```
make TARGET=TestPeriodic1D_MPI.cpp
make run Nproc=4
```

## Two-Dimensional Helmholtz Equation in the Cartesian Coordinate System

We solve the following equation using the Finite Difference Method(FDM)

$$\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} - \lambda u = f,\quad \lambda \geq 0$$

with the periodic boundary condition along the x-direction and following boundary condition along the y-direction
$$\alpha_l \frac{du}{dn}(y_l) + \beta_l u(y_l)=g_l\qquad\text{at the bottom boundary i.e. }y=y_l$$
and
$$\alpha_r \frac{du}{dn}(y_r) + \beta_r u(y_r)=g_r\qquad\text{at the top boundary i.e. }y=y_r$$

Note that $\frac{du}{dn}(y_l) = -\frac{du}{dy}(y_l)$, and $\frac{du}{dn}(y_r) = \frac{du}{dy}(y_r)$.

To build and run the code, invoke the following commands

```
make TARGET=TestNonPeriodic2D_MPI.cpp
make run Nproc=4
```
Here, Nproc denotes the number of processors used to run the code. To use the periodic boundary condition along the y-direction too, build and run the below code

```
make TARGET=TestPeriodic2D_MPI.cpp
make run Nproc=4
```

## Three-Dimensional Helmholtz Equation in the Cartesian Coordinate System

We solve the following equation using the Finite Difference Method(FDM)

$$\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2} - \lambda u = f,\quad \lambda \geq 0$$

with the periodic boundary condition along the x- and y-directions and following boundary condition along the z-direction

$$ \alpha_l \frac{du}{dn}(z_l) + \beta_l u(z_l)=g_l\qquad\text{at the back boundary i.e. } z=z_l $$

and

$$ \alpha_r \frac{du}{dn}(z_r) + \beta_r u(z_r)=g_r\qquad\text{at the front boundary i.e. } z=z_r $$

Note that $\frac{du}{dn}(z_l) = -\frac{du}{dz}(z_l)$, and $\frac{du}{dn}(z_r) = \frac{du}{dz}(z_r)$.

To build and run the code, invoke the following commands

```
make TARGET=TestNonPeriodic3D_MPI.cpp
make run Nproc=4
```
Here, Nproc denotes the number of processors used to run the code. To use the periodic boundary condition along the z-direction too, build and run the below code

```
make TARGET=TestPeriodic3D_MPI.cpp
make run Nproc=4
```
