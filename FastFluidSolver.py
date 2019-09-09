import clr

clr.AddReferenceToFileAndPath(r"C:\Windows\Microsoft.NET\Framework\v4.0.30319\mscorlib.dll")
import System
import System.Collections.Generic
import System.Threading.Tasks
import System.IO
import System.Diagnostics
import System.Text

clr.AddReferenceToFileAndPath(r"C:\Windows\Microsoft.NET\Framework64\v4.0.30319\System.Core.dll")
import System.Linq

clr.AddReferenceToFileAndPath(r"C:\Windows\Microsoft.NET\Framework64\v4.0.30319\WPF\WindowsBase.dll")
import System.Windows

clr.AddReferenceToFileAndPath(r"Z:\ipy2\nuget_dlls\Others\MathNet.Numerics\lib\net461\MathNet.Numerics.dll")
import MathNet.Numerics.LinearAlgebra.Double
import MathNet.Numerics.LinearAlgebra.Storage

# class imports
from System import Math
from System import Array
from System.Collections.Generic import List

class Utilities():
    
    def __init__(self):
        self.EPS = 1*10^(-12)
        
        
    def trilinear_interpolation(self,x=float(),y=float(),z=float(),array=None):
        """ Performs a trilinear interpolation. Method based on code in 
        "Fluid Flow for the Rest of Us: Tutorial of the Marker and Cell Method in Computer Graphics"
        by Cline, Cardon and Egbert.
        
        Parameters:
            "x" - Cell number in x direction
            "y" - Cell number in y direction
            "z" - Cell number in z direction
            "array" - Array to interpolate
        
        Returns:
            Interpolated value
        
        Remarks:
            Cell number here is the cell number on the underlying grid of the array.
            They can be fractional and do not include ghost cells
        """
        imin = Math.Max(Math.Min(int(Math.Floor(x - EPS), array.GetLength(0) - 1)), 0)
        jmin = Math.Max(Math.Min(int(Math.Floor(y - EPS), array.GetLength(1) - 1)), 0)
        kmin = Math.Max(Math.Min(int(Math.Floor(z - EPS), array.GetLength(2) - 1)), 0)
        
        imax = Math.Max(Math.Min(imin + 1, array.GetLength(0) - 1), 0)
        jmax = Math.Max(Math.Min(jmin + 1, array.GetLength(1) - 1), 0)
        kmax = Math.Max(Math.Min(kmin + 1, array.GetLength(2) - 1), 0)
        
        result = (imax - x) * (jmax - y) * (kmax - z) * array[imin, jmin, kmin] + (x - imin) * (jmax - y) * (kmax - z) * array[imax, jmin, kmin] + (imax - x) * (y - jmin) * (kmax - z) * array[imin, jmax, kmin] + (x - imin) * (y - jmin) * (kmax - z) * array[imax, jmax, kmin] + (imax - x) * (jmax - y) * (z - kmin) * array[imin, jmin, kmax] + (x - imin) * (jmax - y) * (z - kmin) * array[imax, jmin, kmax] + (imax - x) * (y - jmin) * (z - kmin) * array[imin, jmax, kmax] + (x - imin) * (y - jmin) * (z - kmin) * array[imax, jmax, kmax]
        
        return result
    
    def compute_L2_difference(x1 = [], x2 = []):
        """
        Computes the pointwise L2 difference between 2 multidimensional arrays.
        Parameters:
            "x1" - Array 1
            "x2" - Array 2
        Returns:
            Normalized pointwise L2 difference between array 1 and array 2.
        """
        diff = 0.0
        self.Sx = x1.GetLength(0) # len(x1)
        self.Sy = x1.GetLength(1) # len(x1)
        self.Sz = x1.GetLength(2) # len(x1)
        
        for i in range(self.Sx):
            for j in range(self.Sy):
                for k in range(self.Sz):
                    diff = diff + Math.Pow(x1[i, j, k] - x2[i, j, k], 2)
        
        diff = Math.Sqrt(diff) / Math.Sqrt(self.Sx * self.Sy * self.Sz)
        return diff
        
        
        def in_domain(coordinate, omega): # coorindate = doubt[], omega = Domain
            """
            Checks if a point is inside the domain.
            Parameters:
                "coordinate" - coordiate (x,y,z) of point
                "omega" - Domain
            Returns:
                True if point is inside fluid domain or on boundary, 
                False if point is outside domain or inside and obstacle/ghost cell
            """
            hx = omega.hx
            hy = omega.hy
            hz = omega.hz
            
            #find domain cell that contains point
            idomain_min = Math.Max(Math.Min(int(Math.Floor((1 - self.EPS) * (coordinate[0] + hx) / hx), omega.Nx - 1)), 0)
            jdomain_min = Math.Max(Math.Min(int(Math.Floor((1 - self.EPS) * (coordinate[1] + hy) / hy), omega.Ny - 1)), 0)
            kdomain_min = Math.Max(Math.Min(int(Math.Floor((1 - self.EPS) * (coordinate[2] + hz) / hz), omega.Nz - 1)), 0)
            
            idomain_max = Math.Max(Math.Min(int(Math.Floor((1 + self.EPS) * (coordinate[0] + hx) / hx), omega.Nx - 1)), 0)
            jdomain_max = Math.Max(Math.Min(int(Math.Floor((1 + self.EPS) * (coordinate[1] + hy) / hy), omega.Ny - 1)), 0)
            kdomain_max = Math.Max(Math.Min(int(Math.Floor((1 + self.EPS) * (coordinate[2] + hz) / hz), omega.Nz - 1)), 0)
            
            possibleCelli = List<int>()
            possibleCellj = List<int>()
            possibleCellk = List<int>()
            
            possibleCelli.Add(idomain_min)
            possibleCellj.Add(jdomain_min)
            possibleCellk.Add(kdomain_min)
            
            if (idomain_min != idomain_max):
                possibleCelli.Add(idomain_max)
                
            if (jdomain_min != jdomain_max):
                possibleCellj.Add(jdomain_max)
            
            if (kdomain_min != kdomain_max):
                possibleCellk.Add(kdomain_max)
            
            
            indomain = False
            
            for i in possibleCelli:
                for j in possibleCellj:
                    for k in possibleCellk:
                        if omega.obstacle_cells[i,j,k] == 0:
                            indomain = true
                            break
            return indomain
        
        def calculate_errors(fs, omega, t, component):
            """
            Computes L2 and L-infinity error between FFD approximation and exact solution.
            
            Parameters:
                "fs">FluidSolver containing FFD solution at time t
                "omega">Domain on which exact solution is given
                "t">Time
                "component">Component to evaluate
            
            Returns (Yields):
                "err_l2">Normalized pointwise L2 error
                "err_inf">Pointwise L-infinity error
            
            Remarks:
                The component to evaluate must be one of:
                    1 for pressure
                    2 for x component of velocity
                    3 for y component of velocity
                    4 for z component of velocity
                    
            """
            
            de = DataExtractor(omega, fs)
            
            nu = fs.nu
            
            self.Nx = omega.Nx - 1
            self.Ny = omega.Ny - 1
            self.Nz = omega.Nz - 1
            
            err_array = [self.Nx, self.Ny, self.Nz]
            comp_interp = [self.Nx, self.Ny, self.Nz]
            zeros = [self.Nx, self.Ny, self.Nz]
            
            for i in range(self.Nx):
                for j in range(self.Ny):
                    for k in range(self.Nz):
                        x = i * omega.hx
                        y = j * omega.hy
                        z = k * omega.hz
                        
                        coordinate = [ x, y, z ]
                        velocity_interp = de.get_velocity(x, y, z)
                        
                        u_exact, v_exact, w_exact, p_exact = omega.exact_solution(coordinate, nu, t) #?
                        
                        if component == 1:
                            comp_interp[i, j, k] = de.get_pressure(x, y, z)
                            err_array[i, j, k] = Math.Abs(de.get_pressure(x, y, z) - p_exact)
                        elif component == 2:
                            comp_interp[i, j, k] = velocity_interp[0]
                            err_array[i, j, k] = Math.Abs(velocity_interp[0] - u_exact)
                        elif component == 3:
                            comp_interp[i, j, k] = velocity_interp[1]
                            err_array[i, j, k] = Math.Abs(velocity_interp[1] - v_exact)
                        elif component == 4:
                            comp_interp[i, j, k] = velocity_interp[2]
                            err_array[i, j, k] = Math.Abs(velocity_interp[2] - w_exact)
                        else:
                            print ("error")
                        
            err_l2 = Utilities.compute_L2_difference(err_array, zeros) # L2 norm of errors
            err_inf = err_array.Cast<double>().Max()
        return (err_12,err_inf) #?

class Domain():
    
    def __init__(self):
        self.Nx = int()
        self.Ny = int()
        self.Nz = int()
        
        self.hx = float()
        self.hy = float()
        self.hz = float()
        
        self.length_x = float()
        self.length_y = float()
        self.length_z = float()
        
        # flag to indicate if cell borders a boundary
        self.boundary_cells = [] #{ get; protected set; } int[,,]
        
        # flag to indicate if cell is part of an obstacle
        self.obstacle_cells = [] #{ get; protected set; } int[,,]
        
        self.boundary_normal_x = [] #{ get; protected set; } int[,,] flag to indicate if boundary at cell is normal to x direction
        self.boundary_normal_y = [] #{ get; protected set; } int[,,] flag to indicate if boundary at cell is normal to x direction
        self.boundary_normal_z = [] #{ get; protected set; } int[,,] flag to indicate if boundary at cell is normal to x direction
        
        self.outflow_boundary_x = [] #{ get; protected set; } int[,,]
        self.outflow_boundary_y = [] #{ get; protected set; } int[,,]
        self.outflow_boundary_z = [] #{ get; protected set; } int[,,]
        
        self.boundary_u = [] #{ get; protected set; } double[,,] x component of velocity at boundary
        self.boundary_v = [] #{ get; protected set; } double[,,] y component of velocity at boundary
        self.boundary_w = [] #{ get; protected set; } double[,,] z component of velocity at boundary
        
        self.obstacle_list = [] #{ get; protected set; } List<int[]>()
        self.normal_x_list = [] #{ get; protected set; } List<int[]>()
        self.normal_y_list = [] #{ get; protected set; } List<int[]>()
        self.normal_z_list = [] #{ get; protected set; } List<int[]>()
        
    def set_ghost_flags(self):
        """
        Labels "ghost cells" outside domain as obstacles and flags cells adjacent to them as boundary cells.
        """
        # x = 0 and x = length_x boundaries
        for j in range(self.Ny):
            for k in range(self.Nz):
                self.obstacle_cells[0, j, k] = 1
                self.obstacle_cells[self.Nx - 1, j, k] = 1
                
                self.obstacle_list.Add([ 0, j, k ])
                self.obstacle_list.Add([ self.Nx - 1, j, k ])
                
                self.boundary_cells[1, j, k] = 1
                self.boundary_cells[self.Nx - 2, j, k] = 1
                
                self.boundary_normal_x[1, j, k] = -1
                self.boundary_normal_x[self.Nx - 2, j, k] = 1
                
                self.normal_x_list.Add([ 1, j, k, -1 ])
                self.normal_x_list.Add([self.Nx - 2, j, k, 1 ])
        # y = 0 and y = length_y boundaries
        for i in range(self.Nx):
            for k in range(self.Nz):
                obstacle_cells[i, 0, k] = 1
                obstacle_cells[i, self.Ny - 1, k] = 1
                
                obstacle_list.Add([ i, 0, k ])
                obstacle_list.Add([ i, self.Ny - 1, k ])
                
                boundary_cells[i, 1, k] = 1
                boundary_cells[i, self.Ny - 2, k] = 1
                
                boundary_normal_y[i, 1, k] = -1
                boundary_normal_y[i, self.Ny - 2, k] = 1
                
                normal_y_list.Add([ i, 1, k, -1 ])
                normal_y_list.Add([ i, self.Ny - 2, k, 1 ])
        # z = 0 and z = length_z boundaries
        for i in range(self.Nx):
            for j in range(self.Ny):
                obstacle_cells[i, j, 0] = 1
                obstacle_cells[i, j, self.Nz - 1] = 1
                
                obstacle_list.Add([ i, j, 0 ])
                obstacle_list.Add([ i, j, self.Nz - 1 ])
                
                boundary_cells[i, j, 1] = 1
                boundary_cells[i, j, self.Nz - 2] = 1
                
                boundary_normal_z[i, j, 1] = -1
                boundary_normal_z[i, j, self.Nz - 2] = 1
                
                normal_z_list.Add([ i, j, 1, -1 ])
                normal_z_list.Add([ i, j, self.Nz - 2, 1 ])
    
    def set_boundary_flags(self):
        """
        Adds boundary flags to cells adjacent to obstacles.
        """
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    if obstacle_sells[i,j,k] == 1:
                        if obstacle_sells[i+1,j,k] == 0:
                            boundary_cells[i+1,j,k] = 1
                            boundary_normal_x[i + 1, j, k] = -1
                            normal_x_list.Add([ i + 1, j, k, -1 ])
                        if obstacle_sells[i-1,j,k] == 0:
                            boundary_cells[i - 1, j, k] = 1
                            boundary_normal_x[i - 1, j, k] = 1
                            normal_x_list.Add([ i - 1, j, k, 1 ])
                            
                            
                        if (obstacle_cells[i, j + 1, k] == 0):
                            boundary_cells[i, j + 1, k] = 1
                            boundary_normal_y[i, j + 1, k] = -1
                            normal_y_list.Add([ i, j + 1, k, -1 ])
                            
                        if (obstacle_cells[i, j - 1, k] == 0):
                            boundary_cells[i, j - 1, k] = 1
                            boundary_normal_y[i, j - 1, k] = 1
                            normal_y_list.Add([ i, j - 1, k, 1 ])
                        
                        if (obstacle_cells[i, j, k + 1] == 0):
                            boundary_cells[i, j, k + 1] = 1
                            boundary_normal_z[i, j, k + 1] = -1
                            normal_z_list.Add([ i, j, k + 1, -1 ])
                            
                        if (obstacle_cells[i, j, k - 1] == 0):
                            boundary_cells[i, j, k - 1] = 1
                            boundary_normal_z[i, j, k - 1] = 1
                            normal_z_list.Add([ i, j, k - 1, 1 ])
                            
                            
    def add_obstacle(self,xmin,xmax,ymin,ymax,zmin,zmax):
        """
        Adds an obstacle to the domain.
        
        Params:
            "xmin" - minimum x coordinate of obstacle
            "xmax" - maximum x coordinate of obstacle
            "ymin" - minimum y coordinate of obstacle
            "ymax" - maximum y coordinate of obstacle
            "zmin" - minimum z coordinate of obstacle
            "zmax" - maximum z coordinate of obstacle
            
        """
        
        i_start = int(Math.Floor(xmin * (Nx - 2) / length_x))
        i_end = int(Math.Floor(xmax * (Nx - 2) / length_x))
        j_start = int(Math.Floor(ymin * (Ny - 2) / length_y))
        j_end = int(Math.Floor(ymax * (Ny - 2) / length_y))
        k_start = int(Math.Floor(zmin * (Nz - 2) / length_z))
        k_end = int(Math.Floor(zmax * (Nz - 2) / length_z))
        
        for i in range(i_start,i_end):
            for j in range(j_start,j_end):
                for k in range(k_start,k_end):
                    obstacle_cells[i + 1, j + 1, k + 1] = 1
                    # obstacle_list.Add(new int[] { i, j, k })
                    
        set_boundary_flags()
        
    def exact_solution(self, coordinate, nu, t): # yields u_exact, v_exact, w_exact, p_exact
        """
        Computes the exact solution at a point if known, the default returns 0 for everything.
        Params:
            "coordinate">coorindate (x,y,z)
            "nu">fluid viscosity
            "t">time
        Returns/Yields:
            "u_exact">exact x component of velocity
            "v_exact">exact y component of velocity
            "w_exact">exact z component of velocity
            "p_exact">exact pressure
            
        """
        
        
        if u_exact and v_exact and w_exact and p_exact:
            return (u_exact,v_exact,w_exact,p_exact)
        else:
            u_exact = 0
            v_exact = 0
            w_exact = 0
            p_exact = 0
            return (u_exact,v_exact,w_exact,p_exact)
        


class DataExtractor():
    print ("")
    # public class DataExtractor
    # {
        # private Domain omega;
        # private FluidSolver fs;
        
        # /// <summary>
        # /// Constructor
        # /// </summary>
        # /// <param name="omega">Domain</param>
        # /// <param name="fs">FFD simulation</param>
        # public DataExtractor(Domain omega, FluidSolver fs)
        # {
            # this.omega = omega;
            # this.fs = fs;
        # }
        
        # /// <summary>
        # /// Calculate the pressure at a point (x,y,z)
        # /// </summary>
        # /// <param name="x">x coordinate</param>
        # /// <param name="y">x coordinate</param>
        # /// <param name="z">x coordinate</param>
        # /// <returns>pressure</returns>
        # public double get_pressure(double x, double y, double z)
        # {
            # double x_scaled = x / omega.hx;
            # double y_scaled = y / omega.hy;
            # double z_scaled = z / omega.hz;
            
            # return Utilities.trilinear_interpolation(x_scaled + 0.5,
                                    # y_scaled + 0.5, z_scaled + 0.5, fs.p);
        # }
        
        # /// <summary>
        # /// Calculate the velocity (u,v,w) at a point (x,y,z)
        # /// </summary>
        # /// <param name="x">x coordinate</param>
        # /// <param name="y">x coordinate</param>
        # /// <param name="z">x coordinate</param>
        # /// <returns>velocity (u,v,w)</returns>
        # public double[] get_velocity(double x, double y, double z)
        # {
            # double[] velocity = new double[3];
            
            # double x_scaled = x / omega.hx;
            # double y_scaled = y / omega.hy;
            # double z_scaled = z / omega.hz;
            
            # velocity[0] = Utilities.trilinear_interpolation(x_scaled, y_scaled + 0.5, z_scaled + 0.5, fs.u);
            # velocity[1] = Utilities.trilinear_interpolation(x_scaled + 0.5, y_scaled, z_scaled + 0.5, fs.v);
            # velocity[2] = Utilities.trilinear_interpolation(x_scaled + 0.5, y_scaled + 0.5, z_scaled, fs.w);
             
            # return velocity;
        # }
    # }


class FluidSolver(): # This is the main class
    print ("")
    # /// <summary>
    # /// Solves the Navier-Stokes equations using the Fast Fluid Dynamics method
    # /// desribed by Stam in the paper "Stable Fluids".
    # ///
    # /// This implementation uses a staggered grid finite difference method to solve the
    # /// spatial equations and backwards Euler in time. Uses a Jacobi iterative method to solve
    # /// the resulting systems. Supports first or second order semi-Lagranian to resolve
    # /// advection term.
    # ///
    # /// List of possible future improvements:
    # /// 1. Parallelize code, in particular Jacobi solver
    # /// 2. Create lists of obstacle and boundary cells to avoid looping over all cells
    # ///     when applying boundary conditions
    # /// </summary>
    # public class FluidSolver
    # {
        # public struct solver_struct
        # {
            # public int max_iter; //maximum number of iterations for Gauss-Seidel solver
            # public int min_iter; //minimum number of iterations
            # public int backtrace_order;
            # public double tol; //maximum relative error for Gauss-Seidel solver
            # public bool verbose;
        # }
        
        # private solver_struct solver_prams;
        
        # public double[,,] u; // x component of velocity
        # public double[,,] v; // y component of velocity
        # public double[,,] w; // z component of velocity
        # public double[,,] p; // pressure
        
        # private double[,,] u_old; //scratch arrays for velocities
        # private double[,,] v_old;
        # private double[,,] w_old;
        # private double[,,] p_old;
        
        # private double dt;  //time step
        # public int Nx { get; private set; } //number of points in each x coordinate direction
        # public int Ny { get; private set; } //number of points in each y coordinate direction
        # public int Nz { get; private set; } //number of points in each z coordinate direction
        
        # private double hx;   //spacing in x coordinate direction
        # private double hy;   //spacing in x coordinate direction
        # private double hz;   //spacing in x coordinate direction
        # public double nu;    //fluid viscosity
        
        # private Domain omega;
        
        # /// <summary>
        # /// Constructor for the fluid solver class.
        # /// </summary>
        # /// <param name="omega">domain omega</param>
        # /// <param name="dt">time step size</param>
        # /// <param name="nu">viscosity</param>
        # /// <param name="u0">initial x component of velocity</param>
        # /// <param name="v0">initial y component of velocity</param>
        # /// <param name="w0">initial z component of velocity</param>
        # /// <param name="solver_prams">structure containing solver options</param>
        # public FluidSolver(Domain omega, double dt, double nu, double[,,] u0, double[,,] v0,
                # double[,,] w0, solver_struct solver_prams)
        # {
            # Nx = omega.Nx;
            # Ny = omega.Ny;
            # Nz = omega.Nz;
            
            # hx = omega.hx;
            # hy = omega.hy;
            # hz = omega.hz;
            
            # this.dt = dt;
            # this.nu = nu;
            
            # this.omega = omega;
            # this.solver_prams = solver_prams;
            
            # p = new double[Nx, Ny, Nz];
            
            # //set up initial pressure guess
            # for (int i = 0; i < Nx; i++)
            # {
                # for (int j = 0; j < Ny; j++)
                # {
                    # for (int k = 0; k < Nz; k++)
                    # {
                        # double x = (i - 0.5) * hx;
                        # p[i, j, k] = -x;
                        # //p[i,j,k] = 0;
                        # //p[i, j, k] = Math.Sin(Math.PI / omega.length_x * x);
                    # }
                # }
            # }
            
            # u = new double[u0.GetLength(0), u0.GetLength(1), u0.GetLength(2)];
            # v = new double[v0.GetLength(0), v0.GetLength(1), v0.GetLength(2)];
            # w = new double[w0.GetLength(0), w0.GetLength(1), w0.GetLength(2)];
            
            # u_old = new double[u.GetLength(0), u.GetLength(1), u.GetLength(2)];
            # v_old = new double[v.GetLength(0), v.GetLength(1), v.GetLength(2)];
            # w_old = new double[w.GetLength(0), w.GetLength(1), w.GetLength(2)];
            # p_old = new double[p.GetLength(0), p.GetLength(1), p.GetLength(2)];
            
            # Array.Copy(u0, 0, u, 0, u0.Length);
            # Array.Copy(v0, 0, v, 0, v0.Length);
            # Array.Copy(w0, 0, w, 0, w0.Length);
            
            # Array.Copy(u, 0, u_old, 0, u.Length);
            # Array.Copy(v, 0, v_old, 0, v.Length);
            # Array.Copy(w, 0, w_old, 0, w.Length);
        # }
        
        # /// <summary>
        # /// Copy constructor
        # /// </summary>
        # public FluidSolver(FluidSolver old)
        # {
            # Nx = old.Nx;
            # Ny = old.Ny;
            # Nz = old.Nz;
            
            # hx = old.hx;
            # hy = old.hy;
            # hz = old.hz;
            
            # dt = old.dt;
            # nu = old.nu;
            
            # u = old.u;
            # v = old.v;
            # w = old.w;
            # p = old.p;
            
            # u_old = old.u_old;
            # v_old = old.v_old;
            # w_old = old.w_old;
            # p_old = old.p_old;
            
            # solver_prams = old.solver_prams;
        # }
        
        # /// <summary>
        # /// Update velocity by adding forcing tem, these can external forces like gravity or
        # /// buoyancy forces for example.
        # /// </summary>
        # /// <param name="f">forces to add</param>
        # /// <param name="x">velocity component array, one of u, v, w</param>
        # private void add_force(double[,,] f, ref double[,,] x)
        # {
            # int Sx = x.GetLength(0);
            # int Sy = x.GetLength(1);
            # int Sz = x.GetLength(2);
            
            # for (int i = 0; i < Sx; i++)
            # {
                # for (int j = 0; j < Sy; j++)
                # {
                    # for (int k = 0; k < Sz; k++)
                    # {
                        # x[i, j, k] += dt * f[i, j, k];
                    # }
                # }
            # }
        # }
        
        # /// <summary>
        # /// Update velocity/concentration by resolving diffusion term. Solves the diffusion
        # /// equation x_t = L(x) using second order finite difference in space and implicit
        # /// Euler in time.
        # /// </summary>
        # /// <param name="x_old">Old state</param>
        # /// <param name="x_new">New state</param>
        # /// <param name="grid_type">Type of grid used</param>
        # /// <remarks>Since we are using a staggered grid, we have to know what kind of
        # /// grid we are using. The following grid numbers are used throughout the program:
        # /// 1: cell centred (pressure, temperature, concentration)
        # /// 2: centre of faces normal to x direction (x component of velocity)
        # /// 3: centre of faces normal to y direction (y component of velocity)
        # /// 4: centre of faces normal to z direction (z component of velocity)</remarks>
        # private void diffuse(double[,,] x_old, ref double[,,] x_new, int grid_type)
        # {
            # double a = 1 + 2 * nu * dt * (Math.Pow(hx, -2) + Math.Pow(hy, -2) + Math.Pow(hz, -2));
            # double[] c = new double[6];
            
            # double[,,] b = new double[x_old.GetLength(0), x_old.GetLength(1), x_old.GetLength(2)];
            # Array.Copy(x_old, 0, b, 0, x_old.Length);
            
            # c[0] = -dt * nu * Math.Pow(hz, -2);
            # c[1] = -dt * nu * Math.Pow(hy, -2);
            # c[2] = -dt * nu * Math.Pow(hx, -2);
            # c[3] = c[2];
            # c[4] = c[1];
            # c[5] = c[0];
            
            # jacobi_solve(a, c, b, x_old, ref x_new, grid_type);
        # }
        
        # /// <summary>
        # /// Projection step. Solves a Poisson equation for the pressure L(p) = div(u_old)
        # /// using finite difference and then updates the velocities u = u_old - grad(p).
        # /// </summary>
        # private void project()
        # {
            # double[,,] div = new double[Nx - 1, Ny - 1, Nz - 1];
            
            # // Calculate div(u_old) using finite differences
            # for (int i = 0; i < Nx; i++)
            # {
                # for (int j = 0; j < Ny; j++)
                # {
                    # for (int k = 0; k < Nz; k++)
                    # {
                        # if (omega.obstacle_cells[i, j, k] == 0)
                        # {
                            # div[i, j, k] = ((u[i, j, k] - u[i - 1, j, k]) / hx +
                                   # (v[i, j, k] - v[i, j - 1, k]) / hy + (w[i, j, k] - w[i, j, k - 1]) / hz) / dt;
                        # }
                    # }
                # }
            # }
            
            # double a = -2 * (Math.Pow(hx, -2) + Math.Pow(hy, -2) + Math.Pow(hz, -2));
            # double[] c = new double[6];
            
            # c[0] = Math.Pow(hz, -2);
            # c[1] = Math.Pow(hy, -2);
            # c[2] = Math.Pow(hx, -2);
            # c[3] = c[2];
            # c[4] = c[1];
            # c[5] = c[0];
            
            # double[,,] p0 = new double[Nx, Ny, Nz]; // Initial guess for pressure
            # Array.Copy(p, p0, p.Length);
            
            # jacobi_solve(a, c, div, p0, ref p, 1);
            
            # double[] coordinate = new double[3];
            
            # // Update velocity by subtracting grad(p)
            # for (int i = 0; i < u.GetLength(0); i++)
            # {
                # for (int j = 0; j < u.GetLength(1); j++)
                # {
                    # for (int k = 0; k < u.GetLength(2); k++)
                    # {
                        # coordinate[0] = i * hx;
                        # coordinate[1] = (j - 0.5) * hy;
                        # coordinate[2] = (k - 0.5) * hz;
                        
                        # if (Utilities.in_domain(coordinate, omega))
                        # {
                            # u[i, j, k] -= dt * (p[i + 1, j, k] - p[i, j, k]) / hx;
                        # }
                    # }
                # }
            # }
            
            # for (int i = 0; i < v.GetLength(0); i++)
            # {
                # for (int j = 0; j < v.GetLength(1); j++)
                # {
                    # for (int k = 0; k < v.GetLength(2); k++)
                    # {
                        # coordinate[0] = (i - 0.5) * hx;
                        # coordinate[1] = j * hy;
                        # coordinate[2] = (k - 0.5) * hz;
                        
                        # if (Utilities.in_domain(coordinate, omega))
                        # {
                            # v[i, j, k] -= dt * (p[i, j + 1, k] - p[i, j, k]) / hy;
                        # }
                    # }
                # }
            # }
            
            # for (int i = 0; i < w.GetLength(0); i++)
            # {
                # for (int j = 0; j < w.GetLength(1); j++)
                # {
                    # for (int k = 0; k < w.GetLength(2); k++)
                    # {
                        # coordinate[0] = (i - 0.5) * hx;
                        # coordinate[1] = (j - 0.5) * hy;
                        # coordinate[2] = k * hz;
                        
                        # if (Utilities.in_domain(coordinate, omega))
                        # {
                            # w[i, j, k] -= dt * (p[i, j, k + 1] - p[i, j, k]) / hz;
                        # }
                    # }
                # }
            # }
            
            # apply_boundary_conditions_list();
        # }
        
        # /// <summary>
        # /// Advection step. Resolve advection term by using a semi-Langrangian backtracer.
        # /// </summary>
        # /// <param name="x">Updated state</param>
        # /// <param name="x0">Original state</param>
        # /// <param name="velx">x component of velocity</param>
        # /// <param name="vely">y component of velocity</param>
        # /// <param name="velz">z component of velocity</param>
        # /// <param name="grid_type">Grid type, as described in diffusion method</param>
        # private void advect(ref double[,,] x, double[,,] x0, double[,,] velx, double[,,] vely,
                    # double[,,] velz, int grid_type)
        # {
            # int Sx = x.GetLength(0);
            # int Sy = x.GetLength(1);
            # int Sz = x.GetLength(2);
            
            # DataExtractor de = new DataExtractor(omega, this);
            
            # // Loop over every node in x
            # for (int i = 1; i < Sx - 1; i++)
            # {
                # for (int j = 1; j < Sy - 1; j++)
                # {
                    # for (int k = 1; k < Sz - 1; k++)
                    # {
                        # double xCoord, yCoord, zCoord;
                        
                        # double[] velocity0, velocity1;
                        # xCoord = yCoord = zCoord = 0;
                        
                        # // Get coordinate of node
                        # switch (grid_type)
                        # {
                            # case 1:
                                # xCoord = (i - 0.5) * hx;
                                # yCoord = (j - 0.5) * hy;
                                # zCoord = (k - 0.5) * hz;
                                
                                # break;
                                
                            # case 2:
                                # xCoord = i * hx;
                                # yCoord = (j - 0.5) * hy;
                                # zCoord = (k - 0.5) * hz;
                                
                                # break;
                                
                            # case 3:
                                # xCoord = (i - 0.5) * hx;
                                # yCoord = j * hy;
                                # zCoord = (k - 0.5) * hz;
                                
                                # break;
                                
                            # case 4:
                                # xCoord = (i - 0.5) * hx;
                                # yCoord = (j - 0.5) * hy;
                                # zCoord = k * hz;
                                
                                # break;
                        # }
                        
                        # double[] coordinate = new double[] { xCoord, yCoord, zCoord };
                        
                        # if (Utilities.in_domain(coordinate, omega))
                        # {
                            # // Find velocity at node
                            # velocity0 = de.get_velocity(xCoord, yCoord, zCoord);
                            # double[] coordBacktraced = new double[3];
                            
                            # switch (solver_prams.backtrace_order)
                            # {
                                # case 1:
                                
                                    # // Perform linear backtrace to find origin of fluid element
                                    # coordBacktraced[0] = xCoord - dt * (velocity0[0]);
                                    # coordBacktraced[1] = yCoord - dt * (velocity0[1]);
                                    # coordBacktraced[2] = zCoord - dt * (velocity0[2]);
                                    
                                    # if (Utilities.in_domain(coordBacktraced, omega))
                                    # {
                                        # // Set velocity at node to be velocity at backtraced coordinate
                                        # x[i, j, k] = Utilities.trilinear_interpolation(i - (dt / hx) * velocity0[0],
                                                    # j - (dt / hy) * velocity0[1], k - (dt / hz) * velocity0[2], x0);
                                    # }
                                    
                                    # break;
                                    
                                # case 2:
                                
                                    # // Perform two step second order backtrace to find origin of fluid element
                                    # coordBacktraced[0] = xCoord - (dt / 2) * (velocity0[0]);
                                    # coordBacktraced[1] = yCoord - (dt / 2) * (velocity0[1]);
                                    # coordBacktraced[2] = zCoord - (dt / 2) * (velocity0[2]);
                                    
                                    # velocity1 = de.get_velocity(coordBacktraced[0], coordBacktraced[1], coordBacktraced[2]);
                                    
                                    # coordBacktraced[0] -= (dt / 2) * velocity1[0];
                                    # coordBacktraced[1] -= (dt / 2) * velocity1[1];
                                    # coordBacktraced[2] -= (dt / 2) * velocity1[2];
                                    
                                    # velocity1 = de.get_velocity(coordBacktraced[0], coordBacktraced[1], coordBacktraced[2]);
                                    
                                    # if (Utilities.in_domain(coordBacktraced, omega))
                                    # {
                                        # // Set velocity at node to be velocity at backtraced coordinate
                                        # x[i, j, k] = Utilities.trilinear_interpolation(
                                                        # i - (dt / (2 * hx)) * (velocity0[0] + velocity1[0]),
                                                        # j - (dt / (2 * hy)) * (velocity0[1] + velocity1[1]),
                                                        # k - (dt / (2 * hz)) * (velocity0[2] + velocity1[2]), x0);
                                    # }
                                    # break;
                            # }
                        # }
                    # }
                # }
            # }
            
            # apply_boundary_conditions_list();
        # }
        
        # /// <summary>
        # /// Applies the boundary conditions.
        # /// </summary>
        # /// <remarks>The staggered grid makes this the longest part of the code.</remarks>
        # private void apply_boundary_conditions()
        # {
            # // loop over all cells
            # for (int i = 1; i < Nx - 1; i++)
            # {
                # for (int j = 1; j < Ny - 1; j++)
                # {
                    # for (int k = 1; k < Nz - 1; k++)
                    # {
                        # if (omega.obstacle_cells[i, j, k] == 0)
                        # {
                            # if (omega.boundary_cells[i, j, k] == 1)
                            # {
                                # /****************************************************************
                                 # * 6 faces, +x, -x, +y, -y, +z, -z
                                 # *
                                 # * For velocity normal to face, simply prescribe value, for other
                                 # * velocities prescribe average of a point inside the domain and
                                 # * a ghost point outside the domain
                                 # ***************************************************************/
                                # if (omega.boundary_normal_x[i, j, k] == -1)//-x face
                                # {
                                    # p[i - 1, j, k] = p[i, j, k];
                                    
                                    # if (omega.outflow_boundary_x[i, j, k] == 1)
                                    # {
                                        # u[i - 1, j, k] = u[i, j, k];
                                        # v[i - 1, j, k] = v[i, j, k];
                                        # w[i - 1, j, k] = w[i, j, k];
                                    # }
                                    # else
                                    # {
                                        # u[i - 1, j, k] = omega.boundary_u[i - 1, j, k];
                                        # v[i - 1, j, k] = 2 * omega.boundary_v[i - 1, j, k] - v[i, j, k];
                                        # w[i - 1, j, k] = 2 * omega.boundary_w[i - 1, j, k] - w[i, j, k];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == 1)//+x face
                                # {
                                    # p[i + 1, j, k] = p[i, j, k];
                                    
                                    # if (omega.outflow_boundary_x[i, j, k] == 1)
                                    # {
                                        # u[i, j, k] = u[i - 1, j, k];
                                        # v[i + 1, j, k] = v[i, j, k];
                                        # w[i + 1, j, k] = w[i, j, k];
                                    # }
                                    # else
                                    # {
                                        # u[i, j, k] = omega.boundary_u[i, j, k];
                                        # v[i + 1, j, k] = 2 * omega.boundary_v[i + 1, j, k] - v[i, j, k];
                                        # w[i + 1, j, k] = 2 * omega.boundary_w[i + 1, j, k] - w[i, j, k];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_y[i, j, k] == -1)//-y face
                                # {
                                    # p[i, j - 1, k] = p[i, j, k];
                                    
                                    # if (omega.outflow_boundary_y[i, j, k] == 1)
                                    # {
                                        # u[i, j - 1, k] = u[i, j, k];
                                        # v[i, j - 1, k] = v[i, j, k];
                                        # w[i, j - 1, k] = w[i, j, k];
                                    # }
                                    # else
                                    # {
                                        # u[i, j - 1, k] = 2 * omega.boundary_u[i, j - 1, k] - u[i, j, k];
                                        # v[i, j - 1, k] = omega.boundary_v[i, j - 1, k];
                                        # w[i, j - 1, k] = 2 * omega.boundary_w[i, j - 1, k] - w[i, j, k];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_y[i, j, k] == 1)//+y face
                                # {
                                    # p[i, j + 1, k] = p[i, j, k];
                                    
                                    # if (omega.outflow_boundary_y[i, j, k] == 1)
                                    # {
                                        # u[i, j + 1, k] = u[i, j, k];
                                        # v[i, j, k] = v[i, j - 1, k];
                                        # w[i, j + 1, k] = w[i, j, k];
                                    # }
                                    # else
                                    # {
                                        # u[i, j + 1, k] = 2 * omega.boundary_u[i, j + 1, k] - u[i, j, k];
                                        # v[i, j, k] = omega.boundary_v[i, j, k];
                                        # w[i, j + 1, k] = 2 * omega.boundary_w[i, j + 1, k] - w[i, j, k];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_z[i, j, k] == -1)//-z face
                                # {
                                    # p[i, j, k - 1] = p[i, j, k];
                                    
                                    # if (omega.outflow_boundary_z[i, j, k] == 1)
                                    # {
                                        # u[i, j, k - 1] = u[i, j, k];
                                        # v[i, j, k - 1] = v[i, j, k];
                                        # w[i, j, k - 1] = w[i, j, k];
                                    # }
                                    # else
                                    # {
                                        # u[i, j, k - 1] = 2 * omega.boundary_u[i, j, k - 1] - u[i, j, k];
                                        # v[i, j, k - 1] = 2 * omega.boundary_v[i, j, k - 1] - v[i, j, k];
                                        # w[i, j, k - 1] = omega.boundary_w[i, j, k - 1];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_z[i, j, k] == 1)//+z face
                                # {
                                    # p[i, j, k + 1] = p[i, j, k];
                                    
                                    # if (omega.outflow_boundary_z[i, j, k] == 1)
                                    # {
                                        # u[i, j, k + 1] = u[i, j, k];
                                        # v[i, j, k + 1] = v[i, j, k];
                                        # w[i, j, k] = w[i, j, k - 1];
                                    # }
                                    # else
                                    # {
                                        # u[i, j, k + 1] = 2 * omega.boundary_u[i, j, k + 1] - u[i, j, k];
                                        # v[i, j, k + 1] = 2 * omega.boundary_v[i, j, k + 1] - v[i, j, k];
                                        # w[i, j, k] = omega.boundary_w[i, j, k];
                                    # }
                                # }
                                
                                # /********************************************************************
                                 # * 12 edges
                                 # *
                                 # * For velocities normal to a face, but on an edge where that hasn't
                                 # * been assigned yet, prescribe velocity. For velocities tangential to
                                 # * the edge, prescribe an average of 4 points around the edge
                                 # *******************************************************************/
                                 
                                # //-x face
                                # if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    # omega.boundary_normal_y[i, j, k] == -1)
                                # {
                                    # p[i - 1, j - 1, k] = p[i, j, k];
                                    
                                    # u[i - 1, j - 1, k] = 2 * omega.boundary_u[i - 1, j - 1, k] -
                                        # omega.boundary_u[i - 1, j, k];
                                        
                                    # v[i - 1, j - 1, k] = 2 * omega.boundary_v[i - 1, j - 1, k] -
                                            # omega.boundary_v[i, j - 1, k];
                                            
                                    # if (omega.outflow_boundary_x[i, j, k] == 0 &&
                                        # omega.outflow_boundary_y[i, j, k] == 1)
                                    # {
                                        # w[i - 1, j - 1, k] = 4 * omega.boundary_w[i - 1, j - 1, k] -
                                               # w[i - 1, j, k] - 2 * w[i, j, k];
                                    # }
                                    # else
                                    # {
                                        # w[i - 1, j - 1, k] = 4 * omega.boundary_w[i - 1, j - 1, k] +
                                            # w[i, j, k] - 2 * omega.boundary_w[i, j - 1, k] -
                                            # 2 * omega.boundary_w[i - 1, j, k];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    # omega.boundary_normal_y[i, j, k] == 1)
                                # {
                                    # p[i - 1, j + 1, k] = p[i, j, k];
                                    
                                    # u[i - 1, j + 1, k] = 2 * omega.boundary_u[i - 1, j + 1, k] -
                                        # omega.boundary_u[i - 1, j, k];
                                        
                                    # v[i - 1, j, k] = 2 * omega.boundary_v[i - 1, j, k] -
                                        # omega.boundary_v[i, j, k];
                                        
                                    # w[i - 1, j + 1, k] = 4 * omega.boundary_w[i - 1, j + 1, k] +
                                        # w[i, j, k] - 2 * omega.boundary_w[i - 1, j, k] -
                                        # 2 * omega.boundary_w[i, j + 1, k];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    # omega.boundary_normal_z[i, j, k] == -1)
                                # {
                                    # p[i - 1, j, k - 1] = p[i, j, k];
                                    
                                    # u[i - 1, j, k - 1] = 2 * omega.boundary_u[i - 1, j, k - 1] -
                                        # omega.boundary_u[i - 1, j, k];
                                        
                                    # w[i - 1, j, k - 1] = 2 * omega.boundary_w[i - 1, j, k - 1] -
                                        # omega.boundary_w[i, j, k - 1];
                                        
                                    # v[i - 1, j, k - 1] = 4 * omega.boundary_v[i - 1, j, k - 1] +
                                        # v[i, j, k] - 2 * omega.boundary_v[i - 1, j, k] -
                                        # 2 * omega.boundary_v[i, j, k - 1];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    # omega.boundary_normal_z[i, j, k] == 1)
                                # {
                                    # p[i - 1, j, k + 1] = p[i, j, k];
                                    
                                    # u[i - 1, j, k + 1] = 2 * omega.boundary_u[i - 1, j, k + 1] -
                                        # omega.boundary_u[i - 1, j, k];
                                        
                                    # v[i - 1, j, k + 1] = 4 * omega.boundary_v[i - 1, j, k + 1] +
                                        # v[i, j, k] - 2 * omega.boundary_v[i - 1, j, k] -
                                        # 2 * omega.boundary_v[i, j, k + 1];
                                        
                                    # w[i - 1, j, k] = 2 * omega.boundary_w[i - 1, j, k] -
                                        # omega.boundary_w[i, j, k];
                                # }
                                
                                # //+x face
                                # if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    # omega.boundary_normal_y[i, j, k] == -1)
                                # {
                                    # p[i + 1, j - 1, k] = p[i, j, k];
                                    
                                    # v[i + 1, j - 1, k] = 2 * omega.boundary_v[i + 1, j - 1, k] -
                                        # omega.boundary_v[i, j - 1, k];
                                        
                                    # w[i + 1, j - 1, k] = 4 * omega.boundary_w[i + 1, j - 1, k] +
                                        # w[i, j, k] - 2 * omega.boundary_w[i + 1, j, k] -
                                        # 2 * omega.boundary_w[i, j - 1, k];
                                        
                                    # if (omega.outflow_boundary_x[i, j, k] == 1 &&
                                        # omega.outflow_boundary_y[i, j, k] == 1)
                                    # {
                                        # u[i, j - 1, k] = u[i - 1, j - 1, k];
                                    # }
                                    # else
                                    # {
                                        # u[i, j - 1, k] = 2 * omega.boundary_u[i, j - 1, k] -
                                            # omega.boundary_u[i, j, k];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    # omega.boundary_normal_y[i, j, k] == 1)
                                # {
                                    # p[i + 1, j + 1, k] = p[i, j, k];
                                    
                                    # v[i + 1, j, k] = 2 * omega.boundary_v[i + 1, j, k] -
                                        # omega.boundary_v[i, j, k];
                                        
                                    # w[i + 1, j + 1, k] = 4 * omega.boundary_w[i + 1, j + 1, k] +
                                        # w[i, j, k] - 2 * omega.boundary_w[i + 1, j, k] -
                                        # 2 * omega.boundary_w[i, j + 1, k];
                                        
                                    # if (omega.outflow_boundary_x[i, j, k] == 1 &&
                                        # omega.outflow_boundary_y[i, j, k] == 1)
                                    # {
                                        # u[i, j + 1, k] = u[i - 1, j + 1, k];
                                    # }
                                    # else
                                    # {
                                        # u[i, j + 1, k] = 2 * omega.boundary_u[i, j + 1, k] -
                                            # omega.boundary_u[i, j, k];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    # omega.boundary_normal_z[i, j, k] == -1)
                                # {
                                    # p[i + 1, j, k - 1] = p[i, j, k];
                                    
                                    # u[i, j, k - 1] = 2 * omega.boundary_u[i, j, k - 1] -
                                        # omega.boundary_u[i, j, k];
                                        
                                    # w[i + 1, j, k - 1] = 2 * omega.boundary_w[i + 1, j, k - 1] -
                                        # omega.boundary_w[i, j, k - 1];
                                        
                                    # v[i + 1, j, k - 1] = 4 * omega.boundary_v[i + 1, j, k - 1] +
                                        # v[i, j, k] - 2 * omega.boundary_v[i + 1, j, k] -
                                        # 2 * omega.boundary_v[i, j, k - 1];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    # omega.boundary_normal_z[i, j, k] == 1)
                                # {
                                    # p[i + 1, j, k + 1] = p[i, j, k];
                                    
                                    # v[i + 1, j, k + 1] = 4 * omega.boundary_v[i + 1, j, k + 1] +
                                        # v[i, j, k] - 2 * omega.boundary_v[i + 1, j, k] -
                                        # 2 * omega.boundary_v[i, j, k + 1];
                                        
                                    # w[i + 1, j, k] = 2 * omega.boundary_w[i + 1, j, k] -
                                        # omega.boundary_w[i, j, k];
                                        
                                    # if (omega.outflow_boundary_x[i, j, k] == 1 &&
                                        # omega.outflow_boundary_z[i, j, k] == 1)
                                    # {
                                        # u[i, j, k + 1] = u[i - 1, j, k + 1];
                                    # }
                                    # else
                                    # {
                                        # u[i, j, k + 1] = 2 * omega.boundary_u[i, j, k + 1] -
                                            # omega.boundary_u[i, j, k];
                                    # }
                                # }
                                
                                # //y,z faces
                                # if (omega.boundary_normal_y[i, j, k] == -1 &&
                                    # omega.boundary_normal_z[i, j, k] == -1)
                                # {
                                    # p[i, j - 1, k - 1] = p[i, j, k];
                                    
                                    # v[i, j - 1, k - 1] = 2 * omega.boundary_v[i, j - 1, k - 1] -
                                        # omega.boundary_v[i, j - 1, k];
                                        
                                    # w[i, j - 1, k - 1] = 2 * omega.boundary_w[i, j - 1, k - 1] -
                                        # omega.boundary_w[i, j, k - 1];
                                        
                                    # if (omega.outflow_boundary_y[i, j, k] == 1 &&
                                        # omega.outflow_boundary_z[i, j, k] == 0)
                                    # {
                                        # u[i, j - 1, k - 1] = 4 * omega.boundary_u[i, j - 1, k - 1] - u[i, j, k - 1]
                                            # - 2 * u[i, j, k];
                                    # }
                                    # else
                                    # {
                                        # u[i, j - 1, k - 1] = 4 * omega.boundary_u[i, j - 1, k - 1] +
                                            # u[i, j, k] - 2 * omega.boundary_u[i, j, k - 1] -
                                            # 2 * omega.boundary_u[i, j - 1, k];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_y[i, j, k] == -1 &&
                                    # omega.boundary_normal_z[i, j, k] == 1)
                                # {
                                    # p[i, j - 1, k + 1] = p[i, j, k];
                                    
                                    # u[i, j - 1, k + 1] = 4 * omega.boundary_u[i, j - 1, k + 1] +
                                        # u[i, j, k] - 2 * omega.boundary_u[i, j, k + 1] -
                                        # 2 * omega.boundary_u[i, j - 1, k];
                                        
                                    # v[i, j - 1, k + 1] = 2 * omega.boundary_v[i, j - 1, k + 1] -
                                        # omega.boundary_v[i, j - 1, k];
                                        
                                    # w[i, j - 1, k] = 2 * omega.boundary_w[i, j - 1, k] -
                                        # omega.boundary_w[i, j, k];
                                # }
                                
                                # if (omega.boundary_normal_y[i, j, k] == 1 &&
                                    # omega.boundary_normal_z[i, j, k] == -1)
                                # {
                                    # p[i, j + 1, k - 1] = p[i, j, k];
                                    
                                    # v[i, j, k - 1] = 2 * omega.boundary_v[i, j, k - 1] -
                                        # omega.boundary_v[i, j, k];
                                        
                                    # w[i, j + 1, k - 1] = 2 * omega.boundary_w[i, j + 1, k - 1] -
                                        # omega.boundary_w[i, j, k - 1];
                                        
                                    # if (omega.outflow_boundary_y[i, j, k] == 1 &&
                                        # omega.outflow_boundary_z[i, j, k] == 0)
                                    # {
                                        # u[i, j + 1, k - 1] = 4 * omega.boundary_u[i, j + 1, k - 1] - u[i, j, k - 1]
                                            # - 2 * u[i, j, k];
                                    # }
                                    # else
                                    # {
                                        # u[i, j + 1, k - 1] = 4 * omega.boundary_u[i, j + 1, k - 1] +
                                            # u[i, j, k] - 2 * omega.boundary_u[i, j, k - 1] -
                                            # 2 * omega.boundary_u[i, j + 1, k];
                                    # }
                                # }
                                
                                # if (omega.boundary_normal_y[i, j, k] == 1 &&
                                    # omega.boundary_normal_z[i, j, k] == 1)
                                # {
                                    # p[i, j + 1, k + 1] = p[i, j, k];
                                    
                                    # u[i, j + 1, k + 1] = 4 * omega.boundary_u[i, j + 1, k + 1] +
                                        # u[i, j, k] - 2 * omega.boundary_u[i, j, k + 1] -
                                        # 2 * omega.boundary_u[i, j + 1, k];
                                        
                                    # v[i, j, k + 1] = 2 * omega.boundary_v[i, j, k + 1] -
                                        # omega.boundary_v[i, j, k];
                                        
                                    # w[i, j + 1, k] = 2 * omega.boundary_w[i, j + 1, k] -
                                        # omega.boundary_w[i, j, k];
                                # }
                                
                                # /*****************************************************************************
                                 # * 8 corners
                                 # *****************************************************************************/
                                 
                                # if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    # omega.boundary_normal_y[i, j, k] == -1 &&
                                    # omega.boundary_normal_z[i, j, k] == -1)
                                # {
                                    # u[i - 1, j - 1, k - 1] = 4 * omega.boundary_u[i - 1, j - 1, k - 1] +
                                        # u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j - 1, k] -
                                        # 2 * omega.boundary_u[i - 1, j, k - 1];
                                        
                                    # v[i - 1, j - 1, k - 1] = 4 * omega.boundary_v[i - 1, j - 1, k - 1] +
                                        # v[i, j - 1, k] - 2 * omega.boundary_v[i - 1, j - 1, k] -
                                        # 2 * omega.boundary_v[i, j - 1, k - 1];
                                        
                                    # w[i - 1, j - 1, k - 1] = 4 * omega.boundary_w[i - 1, j - 1, k - 1] +
                                        # w[i, j, k - 1] - 2 * omega.boundary_w[i - 1, j, k - 1] -
                                        # 2 * omega.boundary_w[i, j - 1, k - 1];
                                        
                                    # p[i - 1, j - 1, k - 1] = p[i, j, k];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    # omega.boundary_normal_y[i, j, k] == 1 &&
                                    # omega.boundary_normal_z[i, j, k] == -1)
                                # {
                                    # u[i - 1, j + 1, k - 1] = 4 * omega.boundary_u[i - 1, j + 1, k - 1] +
                                        # u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j + 1, k] -
                                        # 2 * omega.boundary_u[i - 1, j, k - 1];
                                        
                                    # w[i - 1, j + 1, k - 1] = 4 * omega.boundary_w[i - 1, j + 1, k - 1] +
                                        # w[i, j, k - 1] - 2 * omega.boundary_w[i - 1, j, k - 1] -
                                        # 2 * omega.boundary_w[i, j + 1, k - 1];
                                        
                                    # p[i - 1, j + 1, k - 1] = p[i, j, k];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    # omega.boundary_normal_y[i, j, k] == 1 &&
                                    # omega.boundary_normal_z[i, j, k] == 1)
                                # {
                                    # u[i - 1, j + 1, k + 1] = 4 * omega.boundary_u[i - 1, j + 1, k + 1] +
                                        # u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j + 1, k] -
                                        # 2 * omega.boundary_u[i - 1, j, k + 1];
                                        
                                    # p[i - 1, j + 1, k + 1] = p[i, j, k];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    # omega.boundary_normal_y[i, j, k] == -1 &&
                                    # omega.boundary_normal_z[i, j, k] == 1)
                                # {
                                    # u[i - 1, j - 1, k + 1] = 4 * omega.boundary_u[i - 1, j - 1, k + 1] +
                                        # u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j - 1, k] -
                                        # 2 * omega.boundary_u[i - 1, j, k + 1];
                                        
                                    # v[i - 1, j - 1, k + 1] = 4 * omega.boundary_v[i - 1, j - 1, k + 1] +
                                        # v[i, j - 1, k] - 2 * omega.boundary_v[i - 1, j - 1, k] -
                                        # 2 * omega.boundary_v[i, j - 1, k + 1];
                                        
                                    # p[i - 1, j - 1, k + 1] = p[i, j, k];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    # omega.boundary_normal_y[i, j, k] == -1 &&
                                    # omega.boundary_normal_z[i, j, k] == -1)
                                # {
                                    # v[i + 1, j - 1, k - 1] = 4 * omega.boundary_v[i + 1, j - 1, k - 1] +
                                        # v[i, j - 1, k] - 2 * omega.boundary_v[i + 1, j - 1, k] -
                                        # 2 * omega.boundary_v[i, j - 1, k - 1];
                                        
                                    # w[i + 1, j - 1, k - 1] = 4 * omega.boundary_w[i + 1, j - 1, k - 1] +
                                        # w[i, j, k - 1] - 2 * omega.boundary_w[i + 1, j, k - 1] -
                                        # 2 * omega.boundary_w[i, j - 1, k - 1];
                                        
                                    # p[i + 1, j - 1, k - 1] = p[i, j, k];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    # omega.boundary_normal_y[i, j, k] == 1 &&
                                    # omega.boundary_normal_z[i, j, k] == -1)
                                # {
                                    # w[i + 1, j + 1, k - 1] = 4 * omega.boundary_w[i + 1, j + 1, k - 1] +
                                        # w[i, j, k - 1] - 2 * omega.boundary_w[i + 1, j, k - 1] -
                                        # 2 * omega.boundary_w[i, j + 1, k - 1];
                                        
                                    # p[i + 1, j + 1, k - 1] = p[i, j, k];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    # omega.boundary_normal_y[i, j, k] == 1 &&
                                    # omega.boundary_normal_z[i, j, k] == 1)
                                # {
                                    # p[i + 1, j + 1, k + 1] = p[i, j, k];
                                # }
                                
                                # if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    # omega.boundary_normal_y[i, j, k] == -1 &&
                                    # omega.boundary_normal_z[i, j, k] == 1)
                                # {
                                    # v[i + 1, j - 1, k + 1] = 4 * omega.boundary_v[i + 1, j - 1, k + 1] +
                                        # v[i, j - 1, k] - 2 * omega.boundary_v[i + 1, j - 1, k] -
                                        # 2 * omega.boundary_v[i, j - 1, k + 1];
                                        
                                    # p[i + 1, j - 1, k + 1] = p[i, j, k];
                                # }
                            # }
                        # }
                    # }
                # }
            # }
        # }
        
        # /// <summary>
        # /// Applies the boundary conditions.
        # /// </summary>
        # /// <remarks>The staggered grid makes this the longest part of the code.</remarks>
        # private void apply_boundary_conditions_list()
        # {
            # /****************************************************************
            # * 6 faces, +x, -x, +y, -y, +z, -z
            # *
            # * For velocity normal to face, simply prescribe value, for other
            # * velocities prescribe average of a point inside the domain and
            # * a ghost point outside the domain
            # ***************************************************************/
            
            # foreach (int[] indices in omega.normal_x_list)
            # {
                # int i = indices[0];
                # int j = indices[1];
                # int k = indices[2];
                # int direction = indices[3];
                
                # if (omega.obstacle_cells[i, j, k] != 1)
                # {
                    # if (direction == -1)//-x face
                    # {
                        # p[i - 1, j, k] = p[i, j, k];
                        
                        # if (omega.outflow_boundary_x[i, j, k] == 1)
                        # {
                            # u[i - 1, j, k] = u[i, j, k];
                            # v[i - 1, j, k] = v[i, j, k];
                            # w[i - 1, j, k] = w[i, j, k];
                        # }
                        # else
                        # {
                            # u[i - 1, j, k] = omega.boundary_u[i - 1, j, k];
                            # v[i - 1, j, k] = 2 * omega.boundary_v[i - 1, j, k] - v[i, j, k];
                            # w[i - 1, j, k] = 2 * omega.boundary_w[i - 1, j, k] - w[i, j, k];
                        # }
                    # }
                    
                    # if (direction == 1)//+x face
                    # {
                        # p[i + 1, j, k] = p[i, j, k];
                        
                        # if (omega.outflow_boundary_x[i, j, k] == 1)
                        # {
                            # u[i, j, k] = u[i - 1, j, k];
                            # v[i + 1, j, k] = v[i, j, k];
                            # w[i + 1, j, k] = w[i, j, k];
                        # }
                        # else
                        # {
                            # u[i, j, k] = omega.boundary_u[i, j, k];
                            # v[i + 1, j, k] = 2 * omega.boundary_v[i + 1, j, k] - v[i, j, k];
                            # w[i + 1, j, k] = 2 * omega.boundary_w[i + 1, j, k] - w[i, j, k];
                        # }
                    # }
                # }
            # }
            
            # foreach (int[] indices in omega.normal_y_list)
            # {
                # int i = indices[0];
                # int j = indices[1];
                # int k = indices[2];
                # int direction = indices[3];
                
                # if (omega.obstacle_cells[i, j, k] != 1)
                # {
                    # if (direction == -1)//-x face
                    # {
                        # p[i, j - 1, k] = p[i, j, k];
                        
                        # if (omega.outflow_boundary_y[i, j, k] == 1)
                        # {
                            # u[i, j - 1, k] = u[i, j, k];
                            # v[i, j - 1, k] = v[i, j, k];
                            # w[i, j - 1, k] = w[i, j, k];
                        # }
                        # else
                        # {
                            # u[i, j - 1, k] = 2 * omega.boundary_u[i, j - 1, k] - u[i, j, k];
                            # v[i, j - 1, k] = omega.boundary_v[i, j - 1, k];
                            # w[i, j - 1, k] = 2 * omega.boundary_w[i, j - 1, k] - w[i, j, k];
                        # }
                    # }
                    
                    # if (direction == 1)//+y face
                    # {
                        # p[i, j + 1, k] = p[i, j, k];
                        
                        # if (omega.outflow_boundary_y[i, j, k] == 1)
                        # {
                            # u[i, j + 1, k] = u[i, j, k];
                            # v[i, j, k] = v[i, j - 1, k];
                            # w[i, j + 1, k] = w[i, j, k];
                        # }
                        # else
                        # {
                            # u[i, j + 1, k] = 2 * omega.boundary_u[i, j + 1, k] - u[i, j, k];
                            # v[i, j, k] = omega.boundary_v[i, j, k];
                            # w[i, j + 1, k] = 2 * omega.boundary_w[i, j + 1, k] - w[i, j, k];
                        # }
                    # }
                # }
            # }
            
            # foreach (int[] indices in omega.normal_z_list)
            # {
                # int i = indices[0];
                # int j = indices[1];
                # int k = indices[2];
                # int direction = indices[3];
                
                # if (omega.obstacle_cells[i, j, k] != 1)
                # {
                    # if (direction == -1)//-x face
                    # {
                        # p[i, j, k - 1] = p[i, j, k];
                        
                        # if (omega.outflow_boundary_z[i, j, k] == 1)
                        # {
                            # u[i, j, k - 1] = u[i, j, k];
                            # v[i, j, k - 1] = v[i, j, k];
                            # w[i, j, k - 1] = w[i, j, k];
                        # }
                        # else
                        # {
                            # u[i, j, k - 1] = 2 * omega.boundary_u[i, j, k - 1] - u[i, j, k];
                            # v[i, j, k - 1] = 2 * omega.boundary_v[i, j, k - 1] - v[i, j, k];
                            # w[i, j, k - 1] = omega.boundary_w[i, j, k - 1];
                        # }
                    # }
                    
                    # if (direction == 1)//+z face
                    # {
                        # p[i, j, k + 1] = p[i, j, k];
                        
                        # if (omega.outflow_boundary_z[i, j, k] == 1)
                        # {
                            # u[i, j, k + 1] = u[i, j, k];
                            # v[i, j, k + 1] = v[i, j, k];
                            # w[i, j, k] = w[i, j, k - 1];
                        # }
                        # else
                        # {
                            # u[i, j, k + 1] = 2 * omega.boundary_u[i, j, k + 1] - u[i, j, k];
                            # v[i, j, k + 1] = 2 * omega.boundary_v[i, j, k + 1] - v[i, j, k];
                            # w[i, j, k] = omega.boundary_w[i, j, k];
                        # }
                    # }
                # }
            # }
            
            # //TO DO: ADD IN EDGE AND CORNER CASES
        # }
        
        # /// <summary>
        # /// Calculate mass inflow and outflow
        # /// </summary>
        # /// <param name="m_in">Total mass inflow</param>
        # /// <param name="m_out">Total mass outflow</param>
        # private void calculate_mass_flux(out double m_in, out double m_out)
        # {
            # m_in = 0;
            # m_out = 0;
            
            # // loop over all cells
            # for (int i = 1; i < Nx - 1; i++)
            # {
                # for (int j = 1; j < Ny - 1; j++)
                # {
                    # for (int k = 1; k < Nz - 1; k++)
                    # {
                        # if (omega.obstacle_cells[i, j, k] == 0)
                        # {
                            # if (omega.boundary_normal_x[i, j, k] == -1)//-x face
                            # {
                                # if (u[i - 1, j, k] > 0)
                                # {
                                    # m_in += u[i - 1, j, k] * hy * hz;
                                # }
                                # else
                                # {
                                    # m_out -= u[i - 1, j, k] * hy * hz;
                                # }
                            # }
                            
                            # if (omega.boundary_normal_x[i, j, k] == 1)
                            # {
                                # if (u[i, j, k] > 0)
                                # {
                                    # m_out += u[i, j, k] * hy * hz;
                                # }
                                # else
                                # {
                                    # m_in -= u[i, j, k] * hy * hz;
                                # }
                            # }
                            
                            # if (omega.boundary_normal_y[i, j, k] == -1)
                            # {
                                # if (v[i, j - 1, k] > 0)
                                # {
                                    # m_in += v[i, j - 1, k] * hx * hz;
                                # }
                                # else
                                # {
                                    # m_out -= v[i, j - 1, k] * hx * hz;
                                # }
                            # }
                            
                            # if (omega.boundary_normal_y[i, j, k] == 1)
                            # {
                                # if (v[i, j, k] > 0)
                                # {
                                    # m_out += v[i, j, k] * hx * hz;
                                # }
                                # else
                                # {
                                    # m_in -= v[i, j, k] * hx * hz;
                                # }
                            # }
                            
                            # if (omega.boundary_normal_z[i, j, k] == -1)
                            # {
                                # if (w[i, j - 1, k] > 0)
                                # {
                                    # m_in += w[i, j, k - 1] * hx * hy;
                                # }
                                # else
                                # {
                                    # m_out -= w[i, j, k - 1] * hx * hy;
                                # }
                            # }
                            
                            # if (omega.boundary_normal_z[i, j, k] == 1)
                            # {
                                # if (w[i, j, k] > 0)
                                # {
                                    # m_out += w[i, j, k] * hx * hy;
                                # }
                                # else
                                # {
                                    # m_in -= w[i, j, k] * hx * hy;
                                # }
                            # }
                        # }
                    # }
                # }
            # }
        # }
        
        # /// <summary>
        # /// Apply mass correction to outflow as described by Zuo et al in 2010 paper
        # /// "Improvements in FFD Modeling by Using Different Numerical Schemes"
        # /// </summary>
        # /// <param name="alpha">Ratio used in correction factor</param>
        # /// <remarks>alpha = 1 ensures perfect global mass conservation,
        # /// however paper suggests alpha = 0.7 to avoid instability</remarks>
        # private void apply_mass_correction(double alpha)
        # {
            # double m_in, m_out, correction_factor;
            # calculate_mass_flux(out m_in, out m_out);
            
            # if (Math.Abs(m_out) > 1e-8)
            # {
                # correction_factor = 1 + alpha * ((m_in / m_out) - 1);
            # }
            # else
            # {
                # correction_factor = 1;
            # }
            
            # for (int i = 1; i < Nx - 1; i++)
            # {
                # for (int j = 1; j < Ny - 1; j++)
                # {
                    # for (int k = 1; k < Nz - 1; k++)
                    # {
                        # if (omega.obstacle_cells[i, j, k] == 0)
                        # {
                            # if (omega.boundary_normal_x[i, j, k] == -1)
                            # {
                                # u[i - 1, j, k] *= ((u[i - 1, j, k] < 0) ? correction_factor : 1);
                            # }
                            
                            # if (omega.boundary_normal_x[i, j, k] == 1)
                            # {
                                # u[i, j, k] *= ((u[i, j, k] > 0) ? correction_factor : 1);
                            # }
                            
                            # if (omega.boundary_normal_y[i, j, k] == -1)
                            # {
                                # v[i, j - 1, k] *= ((v[i, j - 1, k] < 0) ? correction_factor : 1);
                            # }
                            
                            # if (omega.boundary_normal_y[i, j, k] == 1)
                            # {
                                # v[i, j, k] *= ((v[i, j, k] > 0) ? correction_factor : 1);
                            # }
                            
                            # if (omega.boundary_normal_z[i, j, k] == -1)
                            # {
                                # w[i, j, k - 1] *= ((w[i, j, k - 1] < 0) ? correction_factor : 1);
                            # }
                            
                            # if (omega.boundary_normal_z[i, j, k] == 1)
                            # {
                                # w[i, j, k] *= ((w[i, j, k] > 0) ? correction_factor : 1);
                            # }
                        # }
                    # }
                # }
            # }
        # }
        
        # /// <summary>
        # /// Perform a single time step. Add forces, diffuse, project, advect, project.
        # /// </summary>
        # /// <param name="f_x">x component of forcing term</param>
        # /// <param name="f_y">y component of forcing term</param>
        # /// <param name="f_z">z component of forcing term</param>
        # public void time_step(double[,,] f_x, double[,,] f_y, double[,,] f_z)
        # {
            # apply_boundary_conditions_list();
            
            # add_force(f_x, ref u);
            # add_force(f_y, ref v);
            # add_force(f_z, ref w);
            
            # Array.Copy(u, 0, u_old, 0, u.Length);
            # Array.Copy(v, 0, v_old, 0, v.Length);
            # Array.Copy(w, 0, w_old, 0, w.Length);
            # Array.Copy(p, 0, p_old, 0, p.Length);
            
            # diffuse(u_old, ref u, 2);
            # diffuse(v_old, ref v, 3);
            # diffuse(w_old, ref w, 4);
            
            # project();
            
            # Array.Copy(u, 0, u_old, 0, u.Length);
            # Array.Copy(v, 0, v_old, 0, v.Length);
            # Array.Copy(w, 0, w_old, 0, w.Length);
            # Array.Copy(p, 0, p_old, 0, p.Length);
            
            # advect(ref u, u_old, u_old, v_old, w_old, 2);
            # advect(ref v, v_old, u_old, v_old, w_old, 3);
            # advect(ref w, w_old, u_old, v_old, w_old, 4);
            
            # // apply_mass_correction(0);
            
            # project();
        # }
        
        # /// <summary>
        # /// Solves the sparse banded system given by the finite difference method applied
        # /// to the Poisson or diffusion equation using the iterative Jacobi method.
        # /// </summary>
        # /// <param name="a">coefficient for diagonal entry</param>
        # /// <param name="c">coefficint array other 6 non-zero entries in each row</param>
        # /// <param name="b">right hand side</param>
        # /// <param name="x0">initial guess</param>
        # /// <param name="x1">solution</param>
        # /// <param name="grid_type">grid type as described in diffusion method</param>
        # /// <remarks>The coefficients for the 6 nonzero entries in each row are given in the
        # /// order x[i,j,k-1], x[i,j-1,k], x[i-1,j,k], x[i+1,j,k], x[i,j+1,k, x[i,j,k+1]</remarks>
        # private void jacobi_solve(double a, double[] c, double[,,] b, double[,,] x0, ref double[,,] x1, int grid_type)
        # {
            # Stopwatch stopWatch = new Stopwatch();
            # stopWatch.Start();
            
            # int Sx = x0.GetLength(0);
            # int Sy = x0.GetLength(1);
            # int Sz = x0.GetLength(2);
            
            # int iter = 0;
            # double res = 2 * solver_prams.tol;
            
            # double[] coordinate = new double[3];
            
            # while (iter < solver_prams.min_iter ||
                # (iter < solver_prams.max_iter && res > solver_prams.tol))
            # {
                # /*if (grid_type == 1)
                # {
                    # x1[3, 3, 3] = 0;
                    # x0[3, 3, 3] = 0;
                # }*/
                
                # apply_boundary_conditions_list();
                
                # for (int k = 1; k < Sz - 1; k++)
                # {
                    # for (int j = 1; j < Sy - 1; j++)
                    # {
                        # for (int i = 1; i < Sx - 1; i++)
                        # {
                            # switch (grid_type)
                            # {
                                # case 1:
                                    # coordinate[0] = (i - 0.5) * hx;
                                    # coordinate[1] = (j - 0.5) * hy;
                                    # coordinate[2] = (k - 0.5) * hz;
                                    # break;
                                    
                                # case 2:
                                    # coordinate[0] = i * hx;
                                    # coordinate[1] = (j - 0.5) * hy;
                                    # coordinate[2] = (k - 0.5) * hz;
                                    # break;
                                    
                                # case 3:
                                    # coordinate[0] = (i - 0.5) * hx;
                                    # coordinate[1] = j * hy;
                                    # coordinate[2] = (k - 0.5) * hz;
                                    # break;
                                    
                                # case 4:
                                    # coordinate[0] = (i - 0.5) * hx;
                                    # coordinate[1] = (j - 0.5) * hy;
                                    # coordinate[2] = k * hz;
                                    # break;
                            # }
                            
                            # if (Utilities.in_domain(coordinate, omega))
                            # {
                                # //if (grid_type != 1 || !( i == 3 && j == 3 && k == 3))
                                # {
                                    # x1[i, j, k] = (b[i, j, k] - (c[0] * x0[i, j, k - 1] +
                                        # c[1] * x0[i, j - 1, k] + c[2] * x0[i - 1, j, k] +
                                        # c[3] * x0[i + 1, j, k] + c[4] * x0[i, j + 1, k] +
                                        # c[5] * x0[i, j, k + 1])) / a;
                                # }
                            # }
                        # }
                    # }
                # }
                
                # res = Utilities.compute_L2_difference(x0, x1);
                # iter++;
                
                # Array.Copy(x1, 0, x0, 0, x1.Length);
            # }
            
            # stopWatch.Stop();
            # TimeSpan ts = stopWatch.Elapsed;
            
            # if (solver_prams.verbose)
            # {
                # Console.WriteLine("Jacobi solver completed with residual of {0} in {1} iterations in {2} seconds",
                    # res, iter, ts.TotalSeconds);
            # }
        # }
    # }

class PostProcessor():
    print ("")
    # /// <summary>
    # /// Class for PostProcessing data calculated by the FFD solver in FluidSolver.cs.
    # /// Exports data to VTK files for further analysis. VTK files
    # /// are used often for visualization of scientific data. These files can be read
    # /// by Visit, Paraview and other applications. A toolkit also exists in C++ for analyzing
    # /// this data. More info on VTK is avaiable here: http://www.vtk.org.
    # /// </summary>
    # public class PostProcessor
    # {
        # FluidSolver fs;
        # Domain omega;
        # DataExtractor de;
        
       # /// <summary>
       # /// Constructor
       # /// </summary>
       # /// <param name="fs">FluidSolver containing solution</param>
       # /// <param name="omega">Domain containing geometry information and 
       # /// exact solutions (if known)</param>
        # public PostProcessor(FluidSolver fs, Domain omega)
        # {
            # this.fs = fs;
            # this.omega = omega;
            
            # de = new DataExtractor(omega, fs);
        # }
        
        # /// <summary>
        # /// Exports data to a VTK file. Data can be cell centred or vertex centred. 
        # /// </summary>
        # /// <param name="fname">file name</param>
        # /// <param name="time">time</param>
        # /// <param name="cell_centred">cell centred flag</param>
        # /// <remarks>The file format guide for VTK is avaliable here:
        # /// http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf </remarks>
        # public void export_data_vtk(String fname, double time, bool cell_centred)
        # {
            # int Nx, Ny, Nz;
            
            # double hx = omega.hx;
            # double hy = omega.hy;
            # double hz = omega.hz;
            
            # double[, ,] u_interp;
            # double[, ,] v_interp;
            # double[, ,] w_interp;
            # double[, ,] p_interp;
            # double[, ,] div_interp;
            
            # double[, ,] err_u;
            # double[, ,] err_v;
            # double[, ,] err_w;
            # double[, ,] err_p;
            
            # interpolate_errors(out err_u, out err_v, out err_w, out err_p, time);
            
            # if (cell_centred)
            # {
                # Nx = omega.Nx - 2;
                # Ny = omega.Ny - 2;
                # Nz = omega.Nz - 2;
                
                # interpolate_to_cell_centre_grid(out u_interp, out v_interp, out w_interp, out p_interp, out div_interp);
            # }
            # else
            # {
            
            
                # Nx = omega.Nx - 1;
                # Ny = omega.Ny - 1;
                # Nz = omega.Nz - 1;
                
                # interpolate_to_vertex_grid(out u_interp, out v_interp, out w_interp, out p_interp);
                # div_interp = new double[1, 1, 1]; //have to assign div_interp to avoid compiling errors
            # }
            
            # using (StreamWriter sw = new StreamWriter(fname))
            # {
                # sw.WriteLine("# vtk DataFile Version 3.0");
                # sw.WriteLine("Fast Fluid Dynamics data\n");
                # sw.WriteLine("ASCII");
                # sw.WriteLine("DATASET RECTILINEAR_GRID");
                # sw.WriteLine("FIELD FieldData 1");
                # sw.WriteLine("TIME 1 1 double");
                # sw.WriteLine("{0}", time);
                # sw.WriteLine("DIMENSIONS {0} {1} {2}", omega.Nx - 1, omega.Ny - 1, omega.Nz - 1);
                # sw.WriteLine("X_COORDINATES {0} double", omega.Nx - 1);
                
                # for (int i = 0; i < omega.Nx - 1; i++)
                # {
                    # sw.WriteLine("{0}", hx * i);
                # }
                # sw.WriteLine("Y_COORDINATES {0} double", omega.Ny - 1);
                # for (int i = 0; i < omega.Ny - 1; i++)
                # {
                    # sw.WriteLine("{0}", hy * i);
                # }
                # sw.WriteLine("Z_COORDINATES {0} double", omega.Nz - 1);
                # for (int i = 0; i < omega.Nz - 1; i++)
                # {
                    # sw.WriteLine("{0}", hz * i);
                # }
                
                # if (cell_centred)
                # {
                    # sw.WriteLine("CELL_DATA {0}", Nx * Ny * Nz);
                # }
                # else
                # {
                    # sw.WriteLine("POINT_DATA {0}", Nx * Ny * Nz);
                # }
                
               # sw.WriteLine("VECTORS velocity double");
               # for (int k = 0; k < Nz; k++)
               # {
                   # for (int j = 0; j < Ny; j++)
                   # {
                       # for (int i = 0; i < Nx; i++)
                       # {
                           # sw.WriteLine("{0} {1} {2}", u_interp[i, j, k], v_interp[i, j, k], w_interp[i, j, k]);
                       # }
                   # }
               # }
                
                # sw.WriteLine("SCALARS pressure double {0}", 1);
                # sw.WriteLine("LOOKUP_TABLE default");
                # for (int k = 0; k < Nz; k++)
                # {
                    # for (int j = 0; j < Ny; j++)
                    # {
                        # for (int i = 0; i < Nx; i++)
                        # {
                            # sw.Write("{0} ", p_interp[i, j, k]);
                        # }
                    # }
                # }
                
                # if (cell_centred)
                # {
                    # sw.WriteLine("SCALARS div_u double {0}", 1);
                    # sw.WriteLine("LOOKUP_TABLE default");
                    # for (int k = 0; k < Nz; k++)
                    # {
                        # for (int j = 0; j < Ny; j++)
                        # {
                            # for (int i = 0; i < Nx; i++)
                            # {
                                # sw.Write("{0} ", div_interp[i, j, k]);
                            # }
                        # }
                    # }
                # }
                # else
                # {
                    # sw.WriteLine("SCALARS err_u double {0}", 1);
                    # sw.WriteLine("LOOKUP_TABLE default");
                    # for (int k = 0; k < Nz; k++)
                    # {
                        # for (int j = 0; j < Ny; j++)
                        # {
                            # for (int i = 0; i < Nx; i++)
                            # {
                                # sw.Write("{0} ", err_u[i, j, k]);
                            # }
                        # }
                    # }
                    
                    # sw.WriteLine("SCALARS err_v double {0}", 1);
                    # sw.WriteLine("LOOKUP_TABLE default");
                    # for (int k = 0; k < Nz; k++)
                    # {
                        # for (int j = 0; j < Ny; j++)
                        # {
                            # for (int i = 0; i < Nx; i++)
                            # {
                                # sw.Write("{0} ", err_v[i, j, k]);
                            # }
                        # }
                    # }
                    
                    # sw.WriteLine("SCALARS err_w double {0}", 1);
                    # sw.WriteLine("LOOKUP_TABLE default");
                    # for (int k = 0; k < Nz; k++)
                    # {
                        # for (int j = 0; j < Ny; j++)
                        # {
                            # for (int i = 0; i < Nx; i++)
                            # {
                                # sw.Write("{0} ", err_w[i, j, k]);
                            # }
                        # }
                    # }
                    
                    # sw.WriteLine("SCALARS err_p double {0}", 1);
                    # sw.WriteLine("LOOKUP_TABLE default");
                    # for (int k = 0; k < Nz; k++)
                    # {
                        # for (int j = 0; j < Ny; j++)
                        # {
                            # for (int i = 0; i < Nx; i++)
                            # {
                                # sw.Write("{0} ", err_p[i, j, k]);
                            # }
                        # }
                    # }
                # }
            # }
        # }
        
        
        # /// <summary>
        # /// Exports geometry, including obstacle cells (excluding ghost cells), boundary cells
        # /// and outflow cells to a VTK file
        # /// </summary>
        # /// <param name="fname">file name</param>
        # /// <param name="time">time</param>
       # public void export_geometry_vtk(String fname, double time)
       # {
           # double hx = omega.length_x / (omega.Nx - 2);
           # double hy = omega.length_y / (omega.Ny - 2);
           # double hz = omega.length_z / (omega.Nz - 2);
           
           # using (StreamWriter sw = new StreamWriter(fname))
           # {
               # sw.WriteLine("# vtk DataFile Version 3.0");
               # sw.WriteLine("Fast Fluid Dynamics geometry\n");
               # sw.WriteLine("ASCII");
               # sw.WriteLine("DATASET RECTILINEAR_GRID");
               # sw.WriteLine("FIELD FieldData 1");
               # sw.WriteLine("TIME 1 1 double");
               # sw.WriteLine("{0}", time);
               # sw.WriteLine("DIMENSIONS {0} {1} {2}", omega.Nx + 1, omega.Ny + 1, omega.Nz + 1);
               # sw.WriteLine("X_COORDINATES {0} double", omega.Nx + 1);
               
               # for (int i = 0; i < omega.Nx + 1; i++)
               # {
                   # sw.WriteLine("{0}", hx * (i - 1));
               # }
               # sw.WriteLine("Y_COORDINATES {0} double", omega.Ny + 1);
               # for (int i = 0; i < omega.Ny + 1; i++)
               # {
                   # sw.WriteLine("{0}", hy * (i - 1));
               # }
               # sw.WriteLine("Z_COORDINATES {0} double", omega.Nz + 1);
               # for (int i = 0; i < omega.Nz + 1; i++)
               # {
                   # sw.WriteLine("{0}", hz * (i - 1));
               # }
               
               # sw.WriteLine("CELL_DATA {0}", omega.Nx * omega.Ny * omega.Nz);
               # sw.WriteLine("SCALARS boundary_cells double {0}", 1);
               # sw.WriteLine("LOOKUP_TABLE default");
               # for (int k = 0; k < omega.Nz; k++)
               # {
                   # for (int j = 0; j < omega.Ny; j++)
                   # {
                       # for (int i = 0; i < omega.Nx; i++)
                       # {
                           # sw.Write("{0} ", (double) omega.boundary_cells[i, j, k]);
                       # }
                   # }
               # }
               # sw.WriteLine("SCALARS obstacle_cells double {0}", 1);
               # sw.WriteLine("LOOKUP_TABLE default");
               # for (int k = 0; k < omega.Nz; k++)
               # {
                   # for (int j = 0; j < omega.Ny ; j++)
                   # {
                       # for (int i = 0; i < omega.Nx ; i++)
                       # {
                           # if (i > 0 && i < omega.Nx - 1
                               # && j > 0 && j < omega.Ny - 1
                               # && k > 0 && k < omega.Nz - 1)
                           # {
                               # sw.Write("{0} ", (double)omega.obstacle_cells[i, j, k]);
                           # }
                           # else
                           # {
                               # sw.Write("{0} ", 0);
                               # }
                       # }
                   # }
               # }
               
               # sw.WriteLine("SCALARS boundary_normal_x double {0}", 1);
               # sw.WriteLine("LOOKUP_TABLE default");
               # for (int k = 0; k < omega.Nz; k++)
               # {
                   # for (int j = 0; j < omega.Ny; j++)
                   # {
                       # for (int i = 0; i < omega.Nx; i++)
                       # {
                           # sw.Write("{0} ", (double) omega.boundary_normal_x[i, j, k]);
                       # }
                   # }
               # }
               
               # sw.WriteLine("SCALARS boundary_normal_y double {0}", 1);
               # sw.WriteLine("LOOKUP_TABLE default");
               # for (int k = 0; k < omega.Nz; k++)
               # {
                   # for (int j = 0; j < omega.Ny; j++)
                   # {
                       # for (int i = 0; i < omega.Nx; i++)
                       # {
                           # sw.Write("{0} ", (double) omega.boundary_normal_y[i, j, k]);
                       # }
                   # }
               # }
               
               # sw.WriteLine("SCALARS boundary_normal_z double {0}", 1);
               # sw.WriteLine("LOOKUP_TABLE default");
               # for (int k = 0; k < omega.Nz; k++)
               # {
                   # for (int j = 0; j < omega.Ny; j++)
                   # {
                       # for (int i = 0; i < omega.Nx; i++)
                       # {
                           # sw.Write("{0} ", (double) omega.boundary_normal_z[i, j, k]);
                       # }
                   # }
               # }
               
               # sw.WriteLine("SCALARS outflow_boundary_x double {0}", 1);
               # sw.WriteLine("LOOKUP_TABLE default");
               # for (int k = 0; k < omega.Nz; k++)
               # {
                   # for (int j = 0; j < omega.Ny; j++)
                   # {
                       # for (int i = 0; i < omega.Nx; i++)
                       # {
                           # sw.Write("{0} ", (double) omega.outflow_boundary_x[i, j, k]);
                       # }
                   # }
               # }
               
               # sw.WriteLine("SCALARS outflow_boundary_y double {0}", 1);
               # sw.WriteLine("LOOKUP_TABLE default");
               # for (int k = 0; k < omega.Nz; k++)
               # {
                   # for (int j = 0; j < omega.Ny; j++)
                   # {
                       # for (int i = 0; i < omega.Nx; i++)
                       # {
                           # sw.Write("{0} ", (double)omega.outflow_boundary_y[i, j, k]);
                       # }
                   # }
               # }
               
               # sw.WriteLine("SCALARS outflow_boundary_z double {0}", 1);
               # sw.WriteLine("LOOKUP_TABLE default");
               # for (int k = 0; k < omega.Nz; k++)
               # {
                   # for (int j = 0; j < omega.Ny; j++)
                   # {
                       # for (int i = 0; i < omega.Nx; i++)
                       # {
                           # sw.Write("{0} ", (double)omega.outflow_boundary_z[i, j, k]);
                       # }
                   # }
               # }
           # }
       # }
       
        
      # /// <summary>
      # /// Exports an uninterpolated pressure field to a VTK file
      # /// </summary>
      # /// <param name="fname">file name</param>
      # /// <param name="time">time</param>
       # public void export_uninterpolated_vtk(String fname, double time)
       # {
           # int Nx = fs.p.GetLength(0);
           # int Ny = fs.p.GetLength(1);
           # int Nz = fs.p.GetLength(2);
           
           # double hx = omega.hx;
           # double hy = omega.hy;
           # double hz = omega.hz;
           
           # using (StreamWriter sw = new StreamWriter(fname))
           # {
               # sw.WriteLine("# vtk DataFile Version 3.0");
               # sw.WriteLine("Fast Fluid Dynamics data\n");
               # sw.WriteLine("ASCII");
               # sw.WriteLine("DATASET RECTILINEAR_GRID");
               # sw.WriteLine("FIELD FieldData 1");
               # sw.WriteLine("TIME 1 1 double");
               # sw.WriteLine("{0}", time);
               # sw.WriteLine("DIMENSIONS {0} {1} {2}", omega.Nx + 1, omega.Ny + 1, omega.Nz + 1);
               # sw.WriteLine("X_COORDINATES {0} double", omega.Nx + 1);
               
               # for (int i = 0; i < omega.Nx + 1; i++)
               # {
                   # sw.WriteLine("{0}", (i - 1) * hx);
               # }
               # sw.WriteLine("Y_COORDINATES {0} double", omega.Ny + 1);
               # for (int i = 0; i < omega.Ny + 1; i++)
               # {
                   # sw.WriteLine("{0}", (i - 1) * hy);
               # }
               # sw.WriteLine("Z_COORDINATES {0} double", omega.Nz + 1);
               # for (int i = 0; i < omega.Nz + 1; i++)
               # {
                   # sw.WriteLine("{0}", (i - 1) * hz);
               # }
               
               # sw.WriteLine("CELL_DATA {0}", Nx * Ny * Nz);
               # sw.WriteLine("SCALARS p_uninterp double {0}", 1);
               # sw.WriteLine("LOOKUP_TABLE default");
               # for (int k = 0; k < Nz; k++)
               # {
                   # for (int j = 0; j < Ny; j++)
                   # {
                       # for (int i = 0; i < Nx; i++)
                       # {
                           # sw.Write("{0} ", fs.p[i, j, k]);
                       # }
                   # }
               # }
           # }
       # }
       
       
        # /// <summary>
        # /// Interpolated velocity components, pressure and divergence field to cell centred grid.
        # /// </summary>
        # /// <param name="u_interp">interpolated x component of velocity</param>
        # /// <param name="v_interp">interpolated y component of velocity</param>
        # /// <param name="w_interp">interpolated z component of velocity</param>
        # /// <param name="p_interp">interpolated pressure</param>
        # /// <param name="div_interp">interpolated divergence</param>
        # private void interpolate_to_cell_centre_grid(out double[, ,] u_interp, out double[, ,] v_interp, 
                    # out double[, ,] w_interp, out double[, ,] p_interp, out double[, ,] div_interp)
        # {
            # int Nx = omega.Nx - 2;
            # int Ny = omega.Ny - 2;
            # int Nz = omega.Nz - 2;
            
            # u_interp = new double[Nx, Ny, Nz];
            # v_interp = new double[Nx, Ny, Nz];
            # w_interp = new double[Nx, Ny, Nz];
            # p_interp = new double[Nx, Ny, Nz];
            # div_interp = new double[Nx, Ny, Nz];
            
            # for (int i = 1; i < Nx + 1; i++)
            # {
                # for (int j = 1; j < Ny + 1; j++)
                # {
                    # for (int k = 1; k < Nz + 1; k++)
                    # {
                        # //if (omega.obstacle_cells[i, j, k] == 0)
                        # {
                            # double[] velocity_interp = de.get_velocity((i - 0.5) * omega.hx, (j - 0.5) * omega.hy, (k - 0.5) * omega.hz);
                            # u_interp[i - 1, j - 1, k - 1] = velocity_interp[0];
                            # v_interp[i - 1, j - 1, k - 1] = velocity_interp[1];
                            # w_interp[i - 1, j - 1, k - 1] = velocity_interp[2];
                            
                            # p_interp[i - 1, j - 1, k - 1] = de.get_pressure((i - 0.5) * omega.hx, (j - 0.5) * omega.hy, (k - 0.5) * omega.hz);
                            
                            # div_interp[i - 1, j - 1, k - 1] = (fs.u[i, j, k] - fs.u[i - 1, j, k]) / omega.hx + (fs.v[i, j, k] - fs.v[i, j - 1, k]) / omega.hy +
                                    # (fs.w[i, j, k] - fs.w[i, j, k - 1]) / omega.hz;
                        # }
                    # }
                # }
            # }
        # }
        
        
        # /// <summary>
        # /// Interpolates velocity components and pressure to a vertex centred grid.
        # /// </summary>
        # /// <param name="u_interp">interpolated x component of velocity</param>
        # /// <param name="v_interp">interpolated y component of velocity</param>
        # /// <param name="w_interp">interpolated z component of velocity</param>
        # /// <param name="p_interp">intepolated pressure</param>
        # private void interpolate_to_vertex_grid(out double[, ,] u_interp, out double[, ,] v_interp,
                    # out double[, ,] w_interp, out double[, ,] p_interp)
        # {
            # int Nx = omega.Nx - 1;
            # int Ny = omega.Ny - 1;
            # int Nz = omega.Nz - 1;
            
            # u_interp = new double[Nx, Ny, Nz];
            # v_interp = new double[Nx, Ny, Nz];
            # w_interp = new double[Nx, Ny, Nz];
            # p_interp = new double[Nx, Ny, Nz];
            
            # for (int i = 0; i < Nx; i++)
            # {
                # for (int j = 0; j < Ny; j++)
                # {
                    # for (int k = 0; k < Nz; k++)
                    # {
                        # double[] velocity_interp = de.get_velocity(i * omega.hx, j * omega.hy, k * omega.hz);
                        # u_interp[i, j, k] = velocity_interp[0];
                        # v_interp[i, j, k] = velocity_interp[1];
                        # w_interp[i, j, k] = velocity_interp[2];
                        
                        # p_interp[i, j, k] = de.get_pressure(i * omega.hx, j * omega.hy, k* omega.hz);
                    # }
                # }
            # }
        # }
        
        # /// <summary>
        # /// Interpolates the errors in u, v, w and p on a Domain on which an exact solution is known.
        # /// </summary>
        # /// <param name="err_u">error in x component of velocity</param>
        # /// <param name="err_v">error in ycomponent of velocity</param>
        # /// <param name="err_w">error in z component of velocity</param>
        # /// <param name="err_p">error in pressure</param>
        # /// <param name="t">time</param>
        # private void interpolate_errors(out double[, ,] err_u, out double[, ,] err_v, 
            # out double[, ,] err_w, out double [, ,] err_p, double t)
        # {
            # int Nx = omega.Nx - 1;
            # int Ny = omega.Ny - 1;
            # int Nz = omega.Nz - 1;
            
            
            # err_u = new double[Nx, Ny, Nz];
            # err_v = new double[Nx, Ny, Nz];
            # err_w = new double[Nx, Ny, Nz];
            # err_p = new double[Nx, Ny, Nz];
            
            # for (int i = 0; i < Nx; i++)
            # {
                # for (int j = 0; j < Ny; j++)
                # {
                    # for (int k = 0; k < Nz; k++)
                    # {
                        # double x = i * omega.hx;
                        # double y = j * omega.hy;
                        # double z = k * omega.hz;
                        
                        # double[] coordinate = new double[] { x, y, z };
                        
                        # double[] velocity_interp = de.get_velocity(x, y, z);
                        
                        # double u_exact, v_exact, w_exact, p_exact;
                        
                        # omega.exact_solution(coordinate, fs.nu, t, out u_exact, out v_exact,
                                # out w_exact, out p_exact);
                                
                        # err_u[i, j, k] = velocity_interp[0] - u_exact;
                        
                        # err_v[i, j, k] = velocity_interp[1] - v_exact;
                        
                        # err_w[i, j, k] = velocity_interp[2] - w_exact;
                        
                        # err_p[i, j, k] = de.get_pressure(x, y, z) - p_exact;
                    # }
                # }
            # }
        # }
    # }


