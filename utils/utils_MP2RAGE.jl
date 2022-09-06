using MRIReco
"""
Documentation for acqDataFromMP2RAGE(b::BrukerFile)
"""
function acqDataFromMP2RAGE(b::BrukerFile)
    # extract parameters
    sx = MRIFiles.acqSize(b)[1]
    sy = MRIFiles.pvmMatrix(b)[2]
    sz = MRIFiles.pvmMatrix(b)[3]
    RareFact = parse.(Int,b["RareFactor"])

    # extract position of rawdata
    GradPhaseVect = parse.(Float32,b["GradPhaseVector"])
    GradPhaseVect=GradPhaseVect*sy/2;
    GradPhaseVect=round.(Int,(GradPhaseVect.-minimum(GradPhaseVect)));

    GradSliceVect = parse.(Float32,b["GradSliceVector"])
    GradSliceVect=GradSliceVect*sz/2;
    GradSliceVect=round.(Int,(GradSliceVect.-minimum(GradSliceVect)))

    #=
    mask = zeros(Int,(sy,sz))
    anim = @animate for i = 1 : 100#length(GradPhaseVect)
        mask[GradPhaseVect[i],GradSliceVect[i]]=1;
        heatmap(mask,color = :greys,aspect_ratio = 1,legend = :none , axis=nothing,lims=(1,sz))
    end
    gif(anim, "anim_fps15.gif", fps = 15)
    =#

    # read data
    dtype = Complex{MRIFiles.acqDataType(b)}
    filename = joinpath(b.path, "fid")

    N = MRIFiles.acqSize(b)
    numChannel = parse.(Int,b["PVM_EncNReceivers"])
    profileLength = Int((ceil(N[1]*numChannel*sizeof(dtype)/1024))*1024/sizeof(dtype)) # number of points + zeros

    RareFact = MRIFiles.acqRareFactor(b)
    numSlices = MRIFiles.acqNumSlices(b)
    numEchos = MRIFiles.acqNumEchos(b)
    numEncSteps2 = length(N) == 3 ? N[3] : 1
    numRep = MRIFiles.acqNumRepetitions(b)

    I = open(filename,"r") do fd
        read!(fd,Array{dtype,5}(undef, profileLength,
                                        RareFact,
                                       numEchos,
                                       div(N[2], RareFact),
                                       numRep))[1:N[1]*numChannel,:,:,:,:]
      end;

    I2 = permutedims(I,(1,2,4,3,5)); # put echos in last dimension
    I2 = reshape(I2,sx,numChannel,N[2],numEchos,numRep)


    objOrd = MRIFiles.acqObjOrder(b)
    objOrd = objOrd.-minimum(objOrd)

    gradMatrix = MRIFiles.acqGradMatrix(b)

    offset1 = MRIFiles.acqReadOffset(b)
    offset2 = MRIFiles.acqPhase1Offset(b)
    offset3 = ndims(b) == 2 ? MRIFiles.acqSliceOffset(b) : MRIFiles.pvmEffPhase2Offset(b)

    profiles = MRIBase.Profile[]
    for nR = 1:numRep
        for nEcho=1:numEchos
            for nEnc = 1:N[2]

            counter = EncodingCounters(kspace_encode_step_1=GradPhaseVect[nEnc],
            kspace_encode_step_2=GradSliceVect[nEnc],
            average=0,
            slice=0,
            contrast=nEcho-1,
            phase=0,
            repetition=nR-1,
            set=0,
            segment=0)

            nSl = 1
            G = gradMatrix[:,:,nSl]
                    read_dir = (G[1,1],G[2,1],G[3,1])
                    phase_dir = (G[1,2],G[2,2],G[3,2])
                    slice_dir = (G[1,3],G[2,3],G[3,3])

            # Not sure if the following is correct...
            pos = offset1[nSl]*G[:,1] +
                    offset2[nSl]*G[:,2] +
                    offset3[nSl]*G[:,3]

            position = (pos[1], pos[2], pos[3])

            head = AcquisitionHeader(number_of_samples=sx, idx=counter,
                                        read_dir=read_dir, phase_dir=phase_dir,
                                        slice_dir=slice_dir, position=position,
                                        center_sample=div(sx,2),
                                        available_channels = numChannel, #TODO
                                        active_channels = numChannel)
            traj = Matrix{Float32}(undef,0,0)
            dat = map(ComplexF32, reshape(I2[:,:,nEnc,nEcho,nR],:,numChannel))
            push!(profiles, MRIBase.Profile(head,traj,dat) )
            end
        end
    end
    params = Dict{String,Any}()
    params["trajectory"] = "Cartesian"
    params["encodedSize"] = [sx;sy;sz]
    F = MRIFiles.acqFov(b)
    params["encodedFOV"] = F
    params["receiverChannels"] = numChannel
    params["H1resonanceFrequency_Hz"] = round(Int, parse(Float64,b["SW"])*1000000)
    params["studyID"] = b["VisuStudyId"]
    #params["studyDescription"] = b["ACQ_scan_name"]
    #params["studyInstanceUID"] =
    params["referringPhysicianName"] = MRIFiles.latin1toutf8(b["ACQ_operator"])

    params["patientName"] = b["VisuSubjectName"]

    params["measurementID"] = parse(Int64,b["VisuExperimentNumber"])
    params["seriesDescription"] = b["ACQ_scan_name"]

    params["institutionName"] = MRIFiles.latin1toutf8(b["ACQ_institution"])
    params["stationName"] = b["ACQ_station"]
    params["systemVendor"] = "Bruker"

    params["TR"] = parse.(Float32,b["RecoveryTime"])
    params["delta_TR"] = MRIFiles.acqRepetitionTime(b)
    params["TE"] = parse.(Float32,b["PVM_EchoTime"])
    params["TI"] = parse.(Float32,b["EffectiveTE"])
    params["flipAngle_deg"] = parse.(Int64,[split(b["ExcPulse1"], ", ")[3];split(b["ExcPulse2"], ", ")[3]])
    params["sequence_type"] = MRIFiles.acqProtocolName(b)
    params["echo_spacing"] = MRIFiles.acqInterEchoTime(b)

    raw = RawAcquisitionData(params, profiles);
    return raw
end

"""
Documentation extractKSpace
"""
function extract3DKSpace(acqData::AcquisitionData)
    nx, ny, nz = acqData.encodingSize[1:3]
    numChan, numSl = MRIFiles.numChannels(acqData), MRIFiles.numSlices(acqData)
    numEcho = length(acqData.traj)
    kdata = zeros(ComplexF64, nx * ny * nz, numEcho,numChan)

    for echo = 1:numEcho
        for coil = 1:numChan
            kdata[acqData.subsampleIndices[numEcho], echo,coil] .= kData(acqData, echo, coil, 1)
        end
    end
    kdata = reshape(kdata, nx, ny, nz, numEcho, numChan)
    return kdata
end

"""
Documentation : Crop the central area
"""
function crop(A::Array{T,4}, s::NTuple{3,Int64}) where {T}
    nx, ny, nz = size(A)
    idx_x = div(nx, 2)-div(s[1], 2)+1:div(nx, 2)-div(s[1], 2)+s[1]
    idx_y = div(ny, 2)-div(s[2], 2)+1:div(ny, 2)-div(s[2], 2)+s[2]
    idx_z = div(nz, 2)-div(s[3], 2)+1:div(nz, 2)-div(s[3], 2)+s[3]
    return A[idx_x, idx_y, idx_z,:]
end

"""
Documentation mp2rage -> MP2RAGE from TI1/TI2
"""
function mp2rage(A)
    mp2 = real((conj(A[:,:,:,1,1]).*A[:,:,:,2,1]) ./ (abs.(A[:,:,:,1,1]).^2 + abs.(A[:,:,:,2,1]).^2 ));
    mp2[isnan.(mp2)] .= 0
    return mp2
end

"""
    param2LUT(b::BrukerFile,T1Range=1:10000)

Compute Lookup Table for MP2RAGE

    In :
        - b::BrukerFile
        - T1Range

    Out :
        lookUpTable (NaN .= 0)
"""
function param2LUT(b::BrukerFile,T1Range=1:10000)

    T1 = T1Range # just renaiming

    n = parse.(Int,b["RareFactor"])
    TR = parse.(Float64,b["PVM_RepetitionTime"])
    MP2RAGE_TR = parse.(Float64,b["RecoveryTime"])

    α₁ = parse.(Float64,split(b["ExcPulse1"],", ")[3])
    α₂ = parse.(Float64,split(b["ExcPulse2"],", ")[3])

    TI1,TI2 = parse.(Float64,b["EffectiveTE"])
    E1=vec(exp.(-TR./T1))
    EA=vec(exp.(-(TI1-(n./2-1).*TR)./T1))
    EB=vec(exp.(-(TI2-TI1-n.*TR)./T1))
    EC=vec(exp.(-(MP2RAGE_TR.-(TI2+(n./2).*TR))./T1))

    eff = 1

    # compute mzss =[(B).*(cosd(α₂).*E1).^n+A].*EC+(1-EC);
    B = ((1 .- EA).*(cosd(α₁).*E1).^n + (1 .- E1).*(1 .- (cosd(α₁).*E1).^n)./(1 .- cosd(α₁).*E1)).*EB+(1 .- EB)
    A = (1 .- E1).*((1 .- (cosd(α₂).*E1).^n)./(1 .- cosd(α₂).*E1))

    mzss_num=((B).*(cosd(α₂).*E1).^n+A).*EC+(1 .-EC)
    mzss_denom=(1 .+ eff.*(cosd(α₁).*cosd(α₂)).^n .* exp.(-MP2RAGE_TR./T1))
    mzss=mzss_num./mzss_denom

    # compute GRE1= sind(α₁).*(A.*(cosd(α₁).*E1).^(n./2-1)+B)
    A=-eff.*mzss.*EA+(1 .- EA)
    B=(1 .- E1).*(1 .- (cosd(α₁).*E1).^(n./2-1))./(1 .- cosd(α₁).*E1)
    GRE1=sind(α₁).*(A.*(cosd(α₁).*E1).^(n./2-1)+B)

    # compute GRE2= sind(α₂).*(A-B)
    A=(mzss-(1 .- EC))./(EC.*(cosd(α₂).*E1).^(n./2))
    B=(1 .- E1).*((cosd(α₂).*E1).^(-n./2) .- 1)./(1 .- cosd(α₂).*E1)

    GRE2=sind(α₂).*(A-B)

    lookUpTable=GRE2.*GRE1./(GRE1.*GRE1+GRE2.*GRE2)
    lookUpTable[isnan.(lookUpTable)].= 0
    return lookUpTable
end

"""
    T1FromMP2RAGE(im_MP2::Array{T, 3},b::BrukerFile;T1range=1:10000)

Compute Lookup Table for MP2RAGE

    In :
        - im_MP2::Array{T,3}
        - b::BrukerFile

    Optionnal :
        - T1 range

    Out :
        - T1map
        - lookUpTable
        - T1range
"""
function T1FromMP2RAGE(im_MP2::Array{Float64, 3},b::BrukerFile;T1Range=1:10000)

    # Generate lookUpTable + cut min
    lookUpTable = param2LUT(b::BrukerFile,T1Range)
    maxVal,maxIdx = findmax(lookUpTable)
    T1Range = T1Range[maxIdx:end]
    lookUpTable = lookUpTable[maxIdx:end]

    T1map = zeros(typeof(im_MP2[1]),size(im_MP2))

    for i = 1:length(T1map)
        idxFirst = findfirst(lookUpTable .<= im_MP2[i])

        if !isnothing(idxFirst)
	  T1map[i] = T1Range[idxFirst]
	end

    end

    return T1map, T1Range,lookUpTable

end
