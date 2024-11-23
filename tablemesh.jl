
# mesh for a single continuous beam
# Rotating beam case from Finite Element Solution of Nonlinear Intrinsic Equations for Curved Composite Beams by Hodges, Shang, and Cesnic
# _______
#   / \   
function mesh_table_beam(;L1 = 31.5, #blade length
    L2 = 6.0, #strut lengths
    Nelem1 = 13,
    Nelem2 = 3,
    angleD_joint = 45.0, # angle of joint (0 is straight, 10 is 10 degrees angle of attack out from center of rotation)
    angleD_strut = 30.0, # angle of strut mounting, 0 is perpendicular to the blade
    strutmount_upper = 0.8,
    strutmount_lower = 0.2)
    # vertical=false)

    angle_joint = angleD_joint/360*2*pi

    # Blade section
    mesh_x1 = collect(LinRange(0.0,L1,Nelem1+1))
    # Insert top strut mount point
    if maximum(L1*strutmount_upper .== mesh_x1) # if we are at exactly an existing node, then offset our mount point
        strutmount_upper += 1e-6
    end
    mesh_x1 = sort([mesh_x1;L1*strutmount_upper])

    # pick out the strut mounting indices
    b2s_topidx = findall(x->isapprox(x,L1*strutmount_upper,atol=1e-5*L1),mesh_x1)[1]

    # Insert bottom strut mount point
    if maximum(L1*strutmount_lower .== mesh_x1) # if we are at exactly an existing node, then offset our mount point
        strutmount_lower += 1e-6
    end
    mesh_x1 = sort([mesh_x1;L1*strutmount_lower])

    # pick out the strut mounting indices
    b2s_botidx = findall(x->isapprox(x,L1*strutmount_lower,atol=1e-5*L1),mesh_x1)[1]

    # Add the y and z components
    mesh_y1 = zero(mesh_x1)
    mesh_z1 = zero(mesh_x1).+L2

    # intra-beam connectivity
    conn1 = zeros(length(mesh_z1)-1,2)
    conn1[:,1] = collect(1:length(mesh_z1)-1)
    conn1[:,2] = collect(2:length(mesh_z1))

    # strut 1
    mesh_zstrut1 = collect(LinRange(0.0,L2,Nelem2+1))
    mesh_xstrut1 = zero(mesh_zstrut1).+strutmount_upper*L1
    mesh_ystrut1 = zero(mesh_zstrut1)

    s2b_topidx = length(mesh_x1)+length(mesh_zstrut1)

    conn2 = zeros(length(mesh_zstrut1)-1,2)
    conn2[:,1] = collect(length(mesh_x1)+1:s2b_topidx-1)
    conn2[:,2] = collect(length(mesh_x1)+2:s2b_topidx)

    # strut 2
    mesh_zstrut2 = collect(LinRange(0.0,L2,Nelem2+1))
    mesh_xstrut2 = zero(mesh_zstrut2).+strutmount_lower*L1
    mesh_ystrut2 = zero(mesh_zstrut2)

    s2b_botidx = length(mesh_x1)+length(mesh_zstrut1)+length(mesh_zstrut2)

    conn3 = zeros(length(mesh_zstrut2)-1,2)
    conn3[:,1] = collect(Nelem1+Nelem2*2+2:s2b_botidx-1)
    conn3[:,2] = collect(Nelem1+Nelem2*2+3:s2b_botidx)

    conn = [conn1;conn2;conn3]

    # if vertical
    mesh_z = [mesh_x1;mesh_xstrut1;mesh_xstrut2].- L1*strutmount_lower
    mesh_y = [mesh_y1;mesh_ystrut1;mesh_ystrut2]
    mesh_x = [mesh_z1;mesh_zstrut1;mesh_zstrut2]
    # else
    # mesh_x = [mesh_x1;mesh_xstrut1;mesh_xstrut2].- L1*strutmount_lower
    # mesh_y = [mesh_y1;mesh_ystrut1;mesh_ystrut2]
    # mesh_z = [mesh_z1;mesh_zstrut1;mesh_zstrut2]
    # end

    numNodes = length(mesh_z)
    nodeNum = collect(LinRange(1,numNodes,numNodes))
    numEl = length(conn[:,1])
    elNum = collect(LinRange(1,numEl,numEl))

    # Define Mesh Types
    # Mesh Type: 0-blade 1-tower 2-strut
    meshtype = zeros(Int,numEl)
    meshtype[1:end] .= 0

    #########################
    # .bld equivalent
    #########################

    meshSeg = zeros(Int,3)

    meshSeg[1] = Nelem1
    meshSeg[2] = Nelem2
    meshSeg[3] = Nelem2

    # Not used for the beam case
    structuralSpanLocNorm = []
    structuralNodeNumbers = []
    structuralElNumbers = []
    # end

    mymesh = OWENSFEA.Mesh(nodeNum,numEl,numNodes,mesh_x,mesh_y,mesh_z,elNum,Int.(conn),meshtype,meshSeg,structuralSpanLocNorm,structuralNodeNumbers,structuralElNumbers)

    ####################################
    ##----------Joint Matrix----------##
    ####################################

    #Connect segments
    jointconn = [s2b_topidx b2s_topidx+1
    s2b_botidx b2s_botidx]

    njoint = length(jointconn[:,1])
    ort = OWENS.calculateElementOrientation(mymesh) #TODO: consider getting rid of ort struct for simplification since it isn't used hardly at all
    # println("start")
    Psi_d_joint = zeros(njoint)
    Theta_d_joint = zeros(njoint)
    for jnt = 1:njoint
        elnum_of_joint = findall(x->x==jointconn[jnt,2],ort.elNum) #gives index of the elNum vector which contains the point index we're after. (the elNum vector is a map between point index and element index)
        if length(elnum_of_joint)==0 #Use the other element associated with the joint
            elnum_of_joint = findall(x->x==jointconn[jnt,2]-1,ort.elNum) #TODO: we get away with this since the elements are increasing and there are no two point objects (and this is only a problem for the top of the tower connecting to the blade tops)
        end
        if length(elnum_of_joint)==0
            elnum_of_joint = findall(x->x==jointconn[jnt,1]-1,ort.elNum)
        end
        Psi_d_joint[jnt] = ort.Psi_d[elnum_of_joint[1]]
        Theta_d_joint[jnt] = ort.Theta_d[elnum_of_joint[1]]
    end
    #Joint Types: (0 = weld(fixed), 1=pinned, 2 = hinge joint with axis about slave node element’s e2 axis, 3 = hinge joint axis about slave node element’s e1 axis, 4 = hinge joint axis about slave node element’s e3 axis)

                 #Joint Number,   Joint Connections, Joint Type, Joint Mass, Not Used, Psi_D, Theta_D
    myjoint = [Float64.(1:1:njoint) jointconn zeros(njoint).+2 zeros(njoint) zeros(njoint) Psi_d_joint Theta_d_joint]

    return mymesh, ort, myjoint
end
