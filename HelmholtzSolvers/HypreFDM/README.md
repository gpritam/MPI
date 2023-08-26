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
make TARGET=TestHelmholtz1D_MPI.cpp
make run Nproc=4
```
Here, Nproc denotes the number of processors used to run the code. To use the periodic boundary condition along the x-direction, assign PeriodicLeftRight = true; in the main() function before building the code.

## Two-Dimensional Helmholtz Equation in the Cartesian Coordinate System

We solve the following equation using the Finite Difference Method(FDM)

$$\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} - \lambda u = f,\quad \lambda \geq 0$$

with boundary conditions 
$$\alpha_l \frac{du}{dn}(x_l) + \beta_l u(x_l)=g_l\qquad\text{at the left boundary i.e. }x=x_l$$
$$\alpha_r \frac{du}{dn}(x_r) + \beta_r u(x_r)=g_r\qquad\text{at the right boundary i.e. }x=x_r$$
$$\alpha_b \frac{du}{dn}(y_l) + \beta_b u(y_l)=g_b\qquad\text{at the bottom boundary i.e. }y=y_l$$
and
$$\alpha_t \frac{du}{dn}(y_r) + \beta_t u(y_r)=g_t\qquad\text{at the top boundary i.e. }y=y_r$$

Note that $\frac{du}{dn}(x_l) = -\frac{du}{dx}(x_l)$, $\frac{du}{dn}(x_r) = \frac{du}{dx}(x_r)$, $\frac{du}{dn}(y_l) = -\frac{du}{dy}(y_l)$, and $\frac{du}{dn}(y_r) = \frac{du}{dy}(y_r)$.

To build and run the code, invoke the following commands

```
make TARGET=TestHelmholtz2D_MPI.cpp
make run Nproc=4
```
Here, Nproc denotes the number of processors used to run the code. To use the periodic boundary conditions along x- and y-directions, assign PeriodicLeftRight = true; and PeriodicBottomTop = true; respectively in the main() function before building the code.

## Three-Dimensional Helmholtz Equation in the Cartesian Coordinate System

We solve the following equation using the Finite Difference Method(FDM)

$$\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2} - \lambda u = f,\quad \lambda \geq 0$$

with boundary conditions 
$$\alpha_l \frac{du}{dn}(x_l) + \beta_l u(x_l)=g_l\qquad\text{at the left boundary i.e. }x=x_l$$
$$\alpha_r \frac{du}{dn}(x_r) + \beta_r u(x_r)=g_r\qquad\text{at the right boundary i.e. }x=x_r$$
$$\alpha_b \frac{du}{dn}(y_l) + \beta_b u(y_l)=g_b\qquad\text{at the bottom boundary i.e. }y=y_l$$
$$\alpha_t \frac{du}{dn}(y_r) + \beta_t u(y_r)=g_t\qquad\text{at the top boundary i.e. }y=y_r$$
$$\alpha_{back} \frac{du}{dn}(z_l) + \beta_{back} u(z_l)=g_{back}\qquad\text{at the back boundary i.e. }z=z_l$$
$$\alpha_f \frac{du}{dn}(z_r) + \beta_f u(z_r)=g_f\qquad\text{at the front boundary i.e. }z=z_r$$

Note that $\frac{du}{dn}(x_l) = -\frac{du}{dx}(x_l)$, $\frac{du}{dn}(x_r) = \frac{du}{dx}(x_r)$, $\frac{du}{dn}(y_l) = -\frac{du}{dy}(y_l)$, $\frac{du}{dn}(y_r) = \frac{du}{dy}(y_r)$, $\frac{du}{dn}(z_l) = -\frac{du}{dz}(z_l)$, and $\frac{du}{dn}(z_r) = \frac{du}{dz}(z_r)$.

To build and run the code, invoke the following commands
```
make TARGET=TestHelmholtz3D_MPI.cpp
make run Nproc=4
```
Here, Nproc denotes the number of processors used to run the code. To use the periodic boundary conditions along x-, y- and z-directions, assign PeriodicLeftRight = true; PeriodicBottomTop = true; and PeriodicBackFront = true; respectively in the main() function before building the code.
