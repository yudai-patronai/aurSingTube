import numpy as np
from scipy.sparse import diags, csr_matrix
from scipy.sparse.linalg import spsolve

class PDESolver:
    """
    A computational model skeleton for solving a system of 2D PDEs with mixed boundary conditions:
    - Newton (Robin) conditions on one boundary
    - von Neumann conditions on other boundaries
    
    The class solves a general elliptic PDE of the form:
    ∇·(D∇u) + b·∇u + cu = f
    """
    
    def __init__(self, Lx, Ly, nx, ny, D=1.0, b=(0,0), c=0.0):
        """
        Initialize the PDE solver.
        
        Parameters:
        - Lx, Ly: Domain dimensions in x and y directions
        - nx, ny: Number of grid points in x and y directions
        - D: Diffusion coefficient (can be scalar or tensor)
        - b: Advection coefficients (bx, by)
        - c: Reaction coefficient
        """
        self.Lx = Lx
        self.Ly = Ly
        self.nx = nx
        self.ny = ny
        self.D = D
        self.b = b
        self.c = c
        
        # Grid spacing
        self.dx = Lx / (nx - 1)
        self.dy = Ly / (ny - 1)
        
        # Create grid
        self.x = np.linspace(0, Lx, nx)
        self.y = np.linspace(0, Ly, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')
        
        # Initialize solution array
        self.u = np.zeros((nx, ny))
        
        # Boundary condition parameters
        self.newton_bc_side = 'right'  # default: Newton BC on right side
        self.alpha = 1.0  # Newton BC coefficient
        self.beta = 0.0   # Newton BC coefficient
        self.g = 0.0      # Newton BC external value
        
    def set_newton_bc(self, side, alpha, beta, g):
        """
        Set Newton (Robin) boundary condition parameters.
        
        Parameters:
        - side: 'left', 'right', 'top', or 'bottom'
        - alpha: Coefficient for u term
        - beta: Coefficient for gradient term
        - g: External value
        """
        self.newton_bc_side = side
        self.alpha = alpha
        self.beta = beta
        self.g = g
        
    def set_von_neumann_bc(self, side, value):
        """
        Set von Neumann boundary condition (normal derivative).
        
        Parameters:
        - side: 'left', 'right', 'top', or 'bottom'
        - value: Value of the normal derivative
        """
        if side == 'left':
            self.left_derivative = value
        elif side == 'right':
            self.right_derivative = value
        elif side == 'top':
            self.top_derivative = value
        elif side == 'bottom':
            self.bottom_derivative = value
            
    def set_dirichlet_bc(self, side, value):
        """
        Set Dirichlet boundary condition (for sides not using Newton BC).
        
        Parameters:
        - side: 'left', 'right', 'top', or 'bottom'
        - value: Boundary value (can be constant or array)
        """
        if side == 'left' and self.newton_bc_side != 'left':
            self.u[0, :] = value
        elif side == 'right' and self.newton_bc_side != 'right':
            self.u[-1, :] = value
        elif side == 'top' and self.newton_bc_side != 'top':
            self.u[:, -1] = value
        elif side == 'bottom' and self.newton_bc_side != 'bottom':
            self.u[:, 0] = value
            
    def construct_matrix(self):
        """
        Construct the finite difference matrix for the PDE system.
        Handles mixed boundary conditions.
        """
        n = self.nx * self.ny
        main_diag = np.ones(n) * (-2/self.dx**2 - 2/self.dy**2 + self.c)
        
        # Diffusion coefficients (simplified for constant D)
        Dx = self.D / self.dx**2
        Dy = self.D / self.dy**2
        
        # Advection coefficients
        bx_term = self.b[0] / (2 * self.dx) if self.b[0] != 0 else 0
        by_term = self.b[1] / (2 * self.dy) if self.b[1] != 0 else 0
        
        # x-direction neighbors
        x_plus = np.ones(n-1) * (Dx - bx_term)
        x_minus = np.ones(n-1) * (Dx + bx_term)
        
        # y-direction neighbors
        y_plus = np.ones(n-self.nx) * (Dy - by_term)
        y_minus = np.ones(n-self.nx) * (Dy + by_term)
        
        # Create the matrix
        diagonals = [main_diag, x_plus, x_minus, y_plus, y_minus]
        offsets = [0, 1, -1, self.nx, -self.nx]
        A = diags(diagonals, offsets, shape=(n, n), format='csr')
        
        # Apply boundary conditions
        A = self.apply_boundary_conditions(A)
        
        return A
    
    def apply_boundary_conditions(self, A):
        """
        Modify the matrix to account for boundary conditions.
        """
        n = self.nx * self.ny
        
        # Apply Newton (Robin) BC on the specified side
        if self.newton_bc_side == 'right':
            for j in range(self.ny):
                i = self.nx - 1  # right boundary
                idx = i + j * self.nx
                # Modify matrix for Newton BC: alpha*u + beta*du/dn = g
                A[idx, idx] = self.alpha - 2*self.beta/self.dx  # center coefficient
                A[idx, idx-1] = 2*self.beta/self.dx  # west neighbor
                # Remove east neighbor (outside domain)
                if idx+1 < n:
                    A[idx, idx+1] = 0
                
        elif self.newton_bc_side == 'left':
            for j in range(self.ny):
                i = 0  # left boundary
                idx = i + j * self.nx
                # Modify matrix for Newton BC
                A[idx, idx] = self.alpha + 2*self.beta/self.dx
                A[idx, idx+1] = -2*self.beta/self.dx
                # Remove west neighbor
                if idx-1 >= 0:
                    A[idx, idx-1] = 0
                
        # Similar modifications for top and bottom Newton BCs would go here
        
        # Apply von Neumann BCs on other boundaries
        # (Implementation would depend on which sides have von Neumann conditions)
        
        return A
    
    def construct_rhs(self, f):
        """
        Construct the right-hand side vector, including boundary terms.
        """
        F = f(self.X, self.Y).flatten()
        
        # Apply Newton BC terms
        if self.newton_bc_side == 'right':
            for j in range(self.ny):
                idx = (self.nx - 1) + j * self.nx
                F[idx] = self.g
        elif self.newton_bc_side == 'left':
            for j in range(self.ny):
                idx = 0 + j * self.nx
                F[idx] = self.g
                
        # Apply von Neumann BC terms would go here
        
        return F
    
    def solve(self, f, max_iter=1000, tol=1e-6):
        """
        Solve the PDE system.
        
        Parameters:
        - f: Source term function f(x,y)
        - max_iter: Maximum number of iterations
        - tol: Tolerance for convergence
        """
        A = self.construct_matrix()
        F = self.construct_rhs(f)
        
        # Solve the linear system
        u_flat = spsolve(A, F)
        self.u = u_flat.reshape((self.nx, self.ny))
        
        return self.u
    
    def plot_solution(self):
        """Visualize the solution using matplotlib."""
        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            
            fig = plt.figure(figsize=(10, 6))
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(self.X, self.Y, self.u, cmap='viridis')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('u(x,y)')
            plt.title('Solution of 2D PDE with Mixed Boundary Conditions')
            plt.show()
        except ImportError:
            print("Matplotlib is required for plotting.")
