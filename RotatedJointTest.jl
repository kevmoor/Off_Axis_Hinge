import OWENSFEA
import OWENS

path,_ = splitdir(@__FILE__)
include("$(path)/tablemesh.jl")

import PyPlot
PyPlot.close("all")
PyPlot.pygui(true)
PyPlot.rc("figure", figsize=(4.5, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=4.0)
PyPlot.rc("legend", frameon=true)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.22, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
# PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

tipload = 2e4

#########################################
### Set up mesh
#########################################
Nelem = 12
mymesh,myort,myjoint =  mesh_table_beam(;L1 = 37.5, #blade length
    L2 = 5.0, #strut lengths
    Nelem1 = Nelem,
    Nelem2 = 3,
    angleD_joint = 45.0, # angle of joint (0 is straight, 10 is 10 degrees angle of attack out from center of rotation)
    angleD_strut = 30.0, # angle of strut mounting, 0 is perpendicular to the blade
    strutmount_upper = 0.9,
    strutmount_lower = 0.1)


hinge_angle = 45.17
myort.Psi_d[1:14] .= hinge_angle
myjoint[1,end-1] = hinge_angle
myjoint[2,end-1] = hinge_angle

PyPlot.figure()
for icon = 1:length(mymesh.conn[:,1])
    idx1 = mymesh.conn[icon,1]
    idx2 = mymesh.conn[icon,2]
    PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
    PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    # sleep(0.1)
end

for ijoint = 1:length(myjoint[:,1])
    idx2 = Int(myjoint[ijoint,2])
    idx1 = Int(myjoint[ijoint,3])
    PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
    sleep(0.1)
end

#########################################
### Set up Sectional Properties
#########################################


NuMad_geom_xlscsv_file = "$path/data/SNL34mGeom.csv"
numadIn_bld = OWENS.readNuMadGeomCSV(NuMad_geom_xlscsv_file)

for (i,airfoil) in enumerate(numadIn_bld.airfoil)
    numadIn_bld.airfoil[i] = "$path/airfoils/$airfoil"
end

NuMad_mat_xlscsv_file = "$path/data/SNL34mMaterials.csv"
plyprops_bld = OWENS.readNuMadMaterialsCSV(NuMad_mat_xlscsv_file)

bld_precompoutput,bld_precompinput,lam_U_bld,lam_L_bld,lam_W_bld = OWENS.getOWENSPreCompOutput(numadIn_bld;plyprops=plyprops_bld)
sectionPropsArray_bld = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,30),numadIn_bld,bld_precompoutput;precompinputs=bld_precompinput)
stiff_bld, mass_bld = OWENS.getSectPropsFromOWENSPreComp(LinRange(0,1,30),numadIn_bld,bld_precompoutput;GX=true)

# That was for all the sections, so we want just the center props
centeridx = round(Int,length(sectionPropsArray_bld)/2)
centersecProp = sectionPropsArray_bld[centeridx]
stiff_bld_center = stiff_bld[centeridx]
mass_bld_center = mass_bld[centeridx]
thickness = maximum(bld_precompinput[round(Int,length(sectionPropsArray_bld)/2)].ynode)*bld_precompinput[round(Int,length(sectionPropsArray_bld)/2)].chord
stiff_bld = fill(stiff_bld_center,mymesh.numEl)
mass_bld = fill(mass_bld_center,mymesh.numEl)
sectionPropsArray = fill(centersecProp,mymesh.numEl)
rotationalEffects = ones(mymesh.numEl)

#store data in element object
myel = OWENSFEA.El(sectionPropsArray,myort.Length,myort.Psi_d,myort.Theta_d,myort.Twist_d,rotationalEffects)

######################################
#### Run OWENSFEA Only
#######################################

# Add the tip force
P1 = round(Int,8)
P2 = round(Int,11/13*Nelem)
nodalinputdata = [P1 "F" 1 0.0
P2 "F" 1 0.0]

# nodalTerms = OWENSFEA.readNodalTerms(data = nodalinputdata)
nodalTerms = OWENSFEA.applyConcentratedTerms(mymesh.numNodes, 6;data = nodalinputdata)

Fdof=[P1*6-6+1,P1*6-6+2]#,P2*6-6+1,P2*6-6+2]
Fexternal=[tipload*cosd(hinge_angle),tipload*sind(hinge_angle)]#,tipload*cosd(hinge_angle),tipload*sind(hinge_angle)]

PyPlot.figure("Point load location")
PyPlot.plot(mymesh.x,mymesh.z,".-")
PyPlot.plot(mymesh.x[P1],mymesh.z[P1],"r.")
# PyPlot.plot(mymesh.x[P2],mymesh.z[P2],"r.")

# Strut/table legs are fixed, Hard coded, #NOTE: don't change the discretization without updating this.
pBC = [
20 1 0
20 2 0
20 3 0
20 4 0
20 5 0
20 6 0
16 1 0
16 2 0
16 3 0
16 4 0
16 5 0
16 6 0
]


feamodel = OWENSFEA.FEAModel(;analysisType = "TNB",
    dataOutputFilename = "none",
    joint = myjoint,
    platformTurbineConnectionNodeNumber = 1,
    pBC,
    # nodalTerms,
    nlOn = false,
    gravityOn = false,
    numNodes = mymesh.numNodes,
    RayleighAlpha = 0.0,
    RayleighBeta = 0.0,
    iterationType = "DI")


displ=zeros(mymesh.numNodes*6)
Omega = 0.0
OmegaStart = 0.0

elStorage = OWENSFEA.initialElementCalculations(feamodel,myel,mymesh)
displ,elStrain,staticAnalysisSuccessful,FReaction = OWENSFEA.staticAnalysis(feamodel,mymesh,myel,displ,Omega,OmegaStart,elStorage;reactionNodeNumber=1,Fdof,Fexternal)

# OWENSFEA.modal(feamodel,mymesh,myel;returnDynMatrices=true)
FReactionHist = [FReaction;;FReaction]'
t = [0,1]

# Extract displacements 
disp_x_OW = [displ[i] for i = 1:6:length(displ)]
disp_y_OW = [displ[i] for i = 2:6:length(displ)]
disp_z_OW = [displ[i] for i = 3:6:length(displ)]
disp_curv1_OW = [displ[i] for i = 4:6:length(displ)]
disp_curv2_OW = [displ[i] for i = 5:6:length(displ)]
disp_curv3_OW = [displ[i] for i = 6:6:length(displ)]

gravity = [0, 0, 0.0]

deformFact = 1.0

# Plot displacements
PyPlot.figure("deflection")
PyPlot.plot(mymesh.x,mymesh.z,"k.",label="Undeformed xz")
PyPlot.plot(mymesh.x,mymesh.y,"k.",label="Undeformed zy")
PyPlot.plot(mymesh.y,mymesh.z,"k.",label="Undeformed yz")
PyPlot.plot((mymesh.x+disp_x_OW*deformFact),(mymesh.z+disp_z_OW*deformFact),".",color=plot_cycle[1],label="OWENSxz")
PyPlot.plot((mymesh.x+disp_x_OW*deformFact),(mymesh.y+disp_y_OW*deformFact),".",color=plot_cycle[2],label="OWENSxy")
PyPlot.plot((mymesh.y+disp_y_OW*deformFact),(mymesh.z+disp_z_OW*deformFact),".",color=plot_cycle[3],label="OWENSyz")
PyPlot.legend()
PyPlot.xlabel("x-position (m)")
PyPlot.ylabel("y-position (m)")
# PyPlot.axis("equal")
# PyPlot.xlim([0,0.5])
# PyPlot.savefig("$(path)/figs/displ$(Nelem*2)elem.pdf",transparent = true)

# Extract Strains
eps_x1_OW = [elStrain[i].epsilon_x[1] for i = 1:length(elStrain)]

eps_y_OW = [elStrain[i].epsilon_y[1] for i = 1:length(elStrain)]
eps_y2_OW = [elStrain[i].epsilon_y[2] for i = 1:length(elStrain)]
eps_y3_OW = [elStrain[i].epsilon_y[3] for i = 1:length(elStrain)]
eps_y4_OW = [elStrain[i].epsilon_y[4] for i = 1:length(elStrain)]
eps_y1_OW = (eps_y_OW.+eps_y2_OW.+eps_y3_OW.+eps_y4_OW).*0.25#0.34785484513745385

eps_z_OW = [elStrain[i].epsilon_z[1] for i = 1:length(elStrain)]
eps_z2_OW = [elStrain[i].epsilon_z[2] for i = 1:length(elStrain)]
eps_z3_OW = [elStrain[i].epsilon_z[3] for i = 1:length(elStrain)]
eps_z4_OW = [elStrain[i].epsilon_z[4] for i = 1:length(elStrain)]
eps_z1_OW = (eps_z_OW.+eps_z2_OW.+eps_z3_OW.+eps_z4_OW).*0.25#0.34785484513745385

kappa_x1_OW = [elStrain[i].kappa_x[1] for i = 1:length(elStrain)]

kappa_y1_OW = [elStrain[i].kappa_y[1] for i = 1:length(elStrain)]

kappa_z1_OW = [elStrain[i].kappa_z[1] for i = 1:length(elStrain)]

zforstrainOW = LinRange(0,length(kappa_z1_OW),length(kappa_z1_OW))#mymesh.z[Int.(mymesh.conn[:,1])]
# Load experimental data

Ealuminum = plyprops_bld.plies[6].e1

# PyPlot.figure()
# PyPlot.plot(zforstrainOW,eps_x1_OW,label="epsilon_x")
# PyPlot.legend()
# # PyPlot.savefig("$(path)/figs/epsilon_x$(Nelem*2)elem.pdf",transparent = true)

# PyPlot.figure()
# PyPlot.plot(zforstrainOW,eps_y1_OW,label="epsilon_y")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(zforstrainOW,eps_z1_OW,label="eps_z")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(zforstrainOW,kappa_x1_OW,label="kappa_x")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(zforstrainOW,kappa_y1_OW,label="kappa_y")
# PyPlot.legend()

# PyPlot.figure()
# PyPlot.plot(zforstrainOW,kappa_z1_OW,label="kappa_z")
# PyPlot.legend()

# Freactions
fx = [FReaction[i] for i = 1:6:length(FReaction)]
fy = [FReaction[i] for i = 2:6:length(FReaction)]
fz = [FReaction[i] for i = 3:6:length(FReaction)]
mx = [FReaction[i] for i = 4:6:length(FReaction)]
my = [FReaction[i] for i = 5:6:length(FReaction)]
mz = [FReaction[i] for i = 6:6:length(FReaction)]
zforfreactionow = LinRange(0,length(mz),length(mz))#mymesh.z[Int.(mymesh.conn[:,1])]

PyPlot.figure()
PyPlot.title("$hinge_angle Degrees")
PyPlot.plot(zforfreactionow[1:14],fx[1:14],".-",color=plot_cycle[1],label="Blade")
PyPlot.plot(zforfreactionow[15:17],fx[15:17],".-",color=plot_cycle[2],label="Strut 1")
PyPlot.plot(zforfreactionow[18:end],fx[18:end],".-",color=plot_cycle[3],label="Strut 2")
PyPlot.xlabel("Element Number")
PyPlot.ylabel("Fx (N)")
PyPlot.legend()
PyPlot.savefig("$(path)/fx_angle$hinge_angle.pdf",transparent = true)

PyPlot.figure()
PyPlot.title("$hinge_angle Degrees")
PyPlot.plot(zforfreactionow[1:14],fy[1:14],".-",color=plot_cycle[1],label="Blade")
PyPlot.plot(zforfreactionow[15:17],fy[15:17],".-",color=plot_cycle[2],label="Strut 1")
PyPlot.plot(zforfreactionow[18:end],fy[18:end],".-",color=plot_cycle[3],label="Strut 2")
PyPlot.xlabel("Element Number")
PyPlot.ylabel("Fy (N)")
PyPlot.legend()
PyPlot.savefig("$(path)/fy_angle$hinge_angle.pdf",transparent = true)

PyPlot.figure()
PyPlot.title("$hinge_angle Degrees")
PyPlot.plot(zforfreactionow[1:14],fz[1:14],".-",color=plot_cycle[1],label="Blade")
PyPlot.plot(zforfreactionow[15:17],fz[15:17],".-",color=plot_cycle[2],label="Strut 1")
PyPlot.plot(zforfreactionow[18:end],fz[18:end],".-",color=plot_cycle[3],label="Strut 2")
PyPlot.xlabel("Element Number")
PyPlot.ylabel("Fz (N)")
PyPlot.legend()
PyPlot.savefig("$(path)/fz_angle$hinge_angle.pdf",transparent = true)

PyPlot.figure()
PyPlot.title("$hinge_angle Degrees")
PyPlot.plot(zforfreactionow[1:14],mx[1:14],".-",color=plot_cycle[1],label="Blade")
PyPlot.plot(zforfreactionow[15:17],mx[15:17],".-",color=plot_cycle[2],label="Strut 1")
PyPlot.plot(zforfreactionow[18:end],mx[18:end],".-",color=plot_cycle[3],label="Strut 2")
PyPlot.xlabel("Element Number")
PyPlot.ylabel("Mx (N-m)")
PyPlot.legend()
PyPlot.savefig("$(path)/mx_angle$hinge_angle.pdf",transparent = true)

PyPlot.figure()
PyPlot.title("$hinge_angle Degrees")
PyPlot.plot(zforfreactionow[1:14],my[1:14],".-",color=plot_cycle[1],label="Blade")
PyPlot.plot(zforfreactionow[15:17],my[15:17],".-",color=plot_cycle[2],label="Strut 1")
PyPlot.plot(zforfreactionow[18:end],my[18:end],".-",color=plot_cycle[3],label="Strut 2")
PyPlot.xlabel("Element Number")
PyPlot.ylabel("My (N-m)")
PyPlot.legend()
PyPlot.savefig("$(path)/my_angle$hinge_angle.pdf",transparent = true)

PyPlot.figure()
PyPlot.title("$hinge_angle Degrees")
PyPlot.plot(zforfreactionow[1:14],mz[1:14],".-",color=plot_cycle[1],label="Blade")
PyPlot.plot(zforfreactionow[15:17],mz[15:17],".-",color=plot_cycle[2],label="Strut 1")
PyPlot.plot(zforfreactionow[18:end],mz[18:end],".-",color=plot_cycle[3],label="Strut 2")
PyPlot.xlabel("Element Number")
PyPlot.ylabel("Mz (N-m)")
PyPlot.legend()
PyPlot.savefig("$(path)/mz_angle$hinge_angle.pdf",transparent = true)


# PyPlot.figure("surface_stress")
# PyPlot.plot(zforstrainOW,kappa_y1_OW.*thickness*Ealuminum,label="Surface Bending Stress_OW")
# PyPlot.legend()
# PyPlot.savefig("$(path)/figs/surfaceBendStress$(Nelem*2).pdf",transparent = true)

println("Creating GXBeam Inputs and Saving the 3D mesh to VTK")
system, assembly, sections = OWENS.owens_to_gx(mymesh,myort,myjoint,sectionPropsArray,stiff_bld, mass_bld;VTKmeshfilename="$path/vtk/table1")

t = [0,1.0]
aziHist = [0.0,0.0]
uHist = [displ displ]'
println("Saving VTK time domain files")
OWENS.OWENSFEA_VTK("$path/vtk/table2",t,uHist,system,assembly,sections;scaling=1,azi=aziHist)


println("Saving VTK time domain files")
userPointNames=["EA","EIyy","EIzz","Fx","Fy","Fz","Mx","My","Mz"]
# userPointData[iname,it,ipt] = Float64

# map el props to points using con
userPointData = zeros(length(userPointNames),length(t),mymesh.numNodes)
EA_points = zeros(mymesh.numNodes)
EIyy_points = zeros(mymesh.numNodes)
EIzz_points = zeros(mymesh.numNodes)

fx_points = zeros(mymesh.numNodes)
fy_points = zeros(mymesh.numNodes)
fz_points = zeros(mymesh.numNodes)
mx_points = zeros(mymesh.numNodes)
my_points = zeros(mymesh.numNodes)
mz_points = zeros(mymesh.numNodes)

# Time-invariant data
for iel = 1:length(myel.props)
    # iel = 1
    nodes = mymesh.conn[iel,:]
    EA_points[Int.(nodes)] = myel.props[iel].EA
    EIyy_points[Int.(nodes)] = myel.props[iel].EIyy
    EIzz_points[Int.(nodes)] = myel.props[iel].EIzz
    fx_points[nodes[2]] = fx[iel]
    fy_points[nodes[2]] = fy[iel]
    fz_points[nodes[2]] = fz[iel]
    mx_points[nodes[2]] = mx[iel]
    my_points[nodes[2]] = my[iel]
    mz_points[nodes[2]] = mz[iel]
end

# fill in the big matrix
for it = 1:length(t)

    userPointData[1,it,:] = EA_points
    userPointData[2,it,:] = EIyy_points
    userPointData[3,it,:] = EIzz_points
    userPointData[4,it,:] = fx_points#uHist[it,1:6:end]
    userPointData[5,it,:] = fy_points#uHist[it,2:6:end]
    userPointData[6,it,:] = fz_points#uHist[it,3:6:end]
    userPointData[7,it,:] = mx_points#uHist[it,4:6:end]
    userPointData[8,it,:] = my_points#uHist[it,5:6:end]
    userPointData[9,it,:] = mz_points#uHist[it,6:6:end]
end

azi=aziHist#./aziHist*1e-6
saveName = "$path/vtk/table1"
OWENS.OWENSFEA_VTK(saveName,t,uHist,system,assembly,sections;scaling=1,azi,userPointNames,userPointData)


