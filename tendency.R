#function that determines the change in concentration of a substrate 
# in soil, as a result of diffusion


## pars require the following parameters
## pars=
# c(Diffusion Coefficient "D" [mm^2 s-1],
# fraction of WFPS at saturation, "porewater" [unitless])
# base mineralization, "basemin" [mg mm-3 s-1]
# microbial uptake rate, "uptake" [s-1]
# microbial death rate, "death" [s-2]
# carbon use efficiency, "cue" [unitless]
#

#additional arguments:
#time: time dimension,scalar [seconds]
#state: a vector containing the state variables in 3 d space, concentration [mg mm-3]
#       is created as array(data="row varies fastest",dim=c(nx,ny,nz) -> state[ix,iy,iz]
#       all state variables are concatenated
#       order of state variables: 1) glucose, 2) microbes
#x: Length dimenions [mm]
#y: Depth dimension  [mm]
#z: Vertical dimesntion [mm]
#dx,dy,dz: grid space, scalar [mm]
#Vfree: available Volume dx*dy*dz*porewater, but also accounts for root Volume

# Output
#dc: tendency in concentration and microbes: dc/dt (array)

dcdt = function(time,state,dx,dy,dz,pars,V_free){
  
  nx = dim(V_free)[1]
  ny = dim(V_free)[2]
  nz = dim(V_free)[3]
  
  length = nx*ny*nz #length of the state vector
  c = array(dim=c(nx,ny,nz),data=state[1:length]) #glucose
  m = array(dim=c(nx,ny,nz),data=state[(length+1):(2*length)])
  flux_x = array(dim=c(nx,ny,nz),data=0)
  flux_y = array(dim=c(nx,ny,nz),data=0)
  flux_z = array(dim=c(nx,ny,nz),data=0)
  

  with(as.list(pars),{
    flux_x[2:(nx-1),,] = D/dx*(c[3:nx,,] - 2*c[2:(nx-1),,] + c[1:(nx-2),,])*dy*dz*porewater
    flux_y[,2:(ny-1),] = D/dy*(c[,3:ny,] - 2*c[,2:(ny-1),] + c[,1:(ny-2),])*dx*dz*porewater
    flux_z[,,2:(nz-1)] = D/dz*(c[,,3:nz] - 2*c[,,2:(nz-1)] + c[,,1:(nz-2)])*dx*dy*porewater
    
    #boundary fluxes mass/time
    flux_x[1,,] = D/dx*(c[2,,] - c[1,,])*dy*dz*porewater
    flux_x[nx,,]= D/dx*(c[(nx-1),,] - c[nx,,])*dy*dz*porewater
    flux_y[,1,] = D/dy*(c[,2,] - c[,1,])*dx*dz*porewater
    flux_y[,ny,]= D/dy*(c[,ny-1,] - c[,ny,])*dx*dz*porewater
    flux_z[,,1] = D/dz*(c[,,2] - c[,,1])*dx*dy*porewater
    flux_z[,,nz]=D/dz*(c[,,nz-1] - c[,,nz])*dx*dy*porewater
    
    
    #biology --------- currently no microbial interaction 
    muptake = uptake * c * m 
    mdeath = m*death
    depoly = Vprime*m/(m+Km)
    
    
    dc = (flux_x + flux_y + flux_z)/V_free - muptake + mdeath + depoly
    dm = muptake*cue - mdeath
    
        
    return(list(c(as.vector(dc),as.vector(dm))))
    
  })
  
}
