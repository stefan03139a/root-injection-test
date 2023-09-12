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
