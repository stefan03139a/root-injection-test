#function that determines the change in concentration of a substrate 
# in soil, as a result of diffusion


## pars require the following parameters
## pars=c(Diffusion Coefficient "D" [mm^2 s-1],
## fraction of WFPS at saturation, "porewater" [unitless])

#additional arguments:
#time: time dimension,scalar [seconds]
#state: a vector containing the state variables in 3 d space, concentration [mg mm-3]
#       is created as array(data="row varies fastest",dim=c(nx,ny,nz) -> state[ix,iy,iz]
#x: Length dimenions [mm]
#y: Depth dimension  [mm]
#z: Vertical dimesntion [mm]
#dx,dy,dz: grid space, scalar [mm]
#Vfree: available Volume dx*dy*dz*porewater, but also accounts for root Volume

# Output
#dc: tendency in concentration: dc/dt (array)

dcdt = function(time,state,x,y,z,dx,dy,dz,pars,Vfree){
  with(as.list(c(state,pars)),{
    
    flux_x[2:(nx-1),,] = D/dx*(c[3:nx,,] - 2*c[2:(nx-1),,] + c[1:(nx-2),,])*dy*dz*porewater
    flux_y[,2:(ny-1),] = D/dy*(c[,3:ny,] - 2*c[,2:(ny-1),] + c[,1:(ny-2),])*dx*xz*porewater
    flux_z[,,2:(ny-1)] = D/dz*(c[,,3:nz] - 2*c[,,2:(nz-1)] + c[,,1:(nz-2)])*dx*dy*porewater
    
    #boundary fluxes mass/time
    flux_x[1,,] = D/dx*(c[2,,] - c[1,,])*dy*dz*porewater
    flux_x[nx,,]= D/dx*(c[(nx-1),,] - c[nx,,])*dy*dz*porewater
    flux_y[,1,] = D/dy*(c[,2,] - c[,2,])*dx*dz*porewater
    flux_y[,ny,]= D/dy*(c[,ny-1,] - c[,ny,])*dx*dz*porewater
    flux_z[,,1] = D/dz*(c[,,2] - c[,,1])*dx*dy*porewater
    flux_z[,,nz]=D/dz*(c[,,nz-1] - c[,,nz])*dx*dy*porewater
    
    dc = (flux_x + flux_y + flux_z)/V_free
    
    #biology --------- currently no microbial interaction 
    return(list(dc))
    
  })
  
}


#easiest to make grid sized as close to the diameter of the roo. big as the the diameter of the root:

# the root is attached to the x-pane of the glass
# y is horizontal depth
# z is vertical depth

r0 = 0.5 #enter radius of root

L  =  10#length of box
W  =  10#horizontal depth
H  =  10 #vertical depth
  
#spatial grid setup in a way that the box size is the smallest possible size to 
# the root cross section, and that the root is completely in one box

nx = L %/% (2*r0) # make sure that a single cell is larger than the root dimameter
nx = ifelse((nx %% 2) == 0, nx-1,nx) #make sure it is odd, so that the root in the middle is completely in 1 box
dx = L/nx
x  = seq(0.5*dx,L-0.5*dx,dx)

#same for y direction:
ny = W %/% (2*r0)
dy = W/ny
y  = seq(0.5*dy, W-0.5*dy,dy)

#z direction
nz = round(H/(2*r0))# no constraints needed, just make sure grid size is approximately the sam
dz = H/nz
z = seq(0.5*dz, H-0.5*dz,dz)


Lroot = 4  #root length change number
H0    =  1 #upper depth at which root is placed
H1    = Lroot + H0 #lower depth

porewater = 0.5

#addition
amount_added = 1  #specify one time addition

#calculate where it is added
weight = numeric(nz)
mid    =  (nx+1)/2   #index in x direction where root is
top    = seq(0,nz-1,1)*dz #indicate the top and bottom of z direction
bottom = seq(1,nz,1)*dz

#
idz = which((top>=H0) & (bottom<=H1),arr.ind=TRUE) #the cells that completely cover the root
weight[idz] = 1

idzH0 = idz[1]-1   #cell that receives only part
weight[idzH0] = (bottom[idzH0]-H0)/dz #

idzH1 = idz[length(idz)] + 1 #last cell that only gets part
weight[idzH1]  = (H1 - top[idzH1])/dz

V_straw = r0^2*pi*weight*dz    #volume of root/straw in each box, just a vector
V_free  = array(data=porewater*dx*dy*dz,dim=c(nx,ny,nz))
V_free[mid,1,] = (dx*dy*dz - V_straw)*porewater #special volume in cells where the root is




#parameters:
pars=c(D=1e-6,porewater=porewater)

#initial conditions
c0 = array(data=0,dim=c(nx,ny,nz))
#input from root exudates
input = 0
flux_in = weight*input*dz/L  #mol s-1 in each cell
c0[mid,1,] = flux_in/V_free[mid,1,]
c0[mid,1,] = 1 #right now just placeholder
times=seq(0,3600,600)

a = dcdt(0,c0,x,y,z,dx,dy,dz,pars,Vfree)
stop()

result=ode.3D(y=c0,times=times,nspec=1,dimens=c(nx,ny,nz),method='rk4',x=x,y=y,z=z,dx=dx,dy=dy,dz=dz)


