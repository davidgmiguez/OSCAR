import Pkg;
Pkg.add("CSV")
Pkg.add("Statistics")
Pkg.add("LinearAlgebra")
Pkg.add("Plots")
Pkg.add("DelimitedFiles")
Pkg.add("DataFrames")
Pkg.add("StatsKit")
Pkg.add("Polynomials")
Pkg.add("Tables")
# using Polynomials
using StatsKit
using CSV
using Statistics
using LinearAlgebra
using Plots
using DelimitedFiles
using DataFrames
using Tables

function data2d_chunking(file, xind, yind, areaind, zind, angind, aind, bind, typeind, nCh)
    # data = CSV.read(file, delim='\t', copycols=true)
    data = CSV.read(file, delim='\t', DataFrame);
    data2d = Vector{Array{Float64}}()
    numberObj2D = Vector{}();
    Zi = zeros(size(data)[1],11+nCh); lookcount = 0; zt = data[1,zind]; indexcount = 1;
    for i in 1:size(data)[1]
        # println("i = ", i)
        if ((data[i,zind]) > zt)
            Zi = Zi[1:lookcount, :];
            push!(data2d,Zi)
            push!(numberObj2D, lookcount);
            # println("rest lookcount")
            lookcount = 0; indexcount = indexcount + 1;
            zt = data[i,zind]
            Zi = zeros(size(data)[1],11+nCh);
        end
        lookcount += 1
        Zi[lookcount,1] = data[i,xind] #Xc
        Zi[lookcount,2] = data[i,yind] #Yc
        Zi[lookcount,3] = data[i,areaind] #Area
        Zi[lookcount,4] = data[i,angind] #Angle
        Zi[lookcount,5] = data[i,zind] #Z
        Zi[lookcount,6] = data[i,aind] #Major exe
        Zi[lookcount,7] = data[i,bind] #Minor exe
        Zi[lookcount,8] = 0 #Tag of object
        Zi[lookcount,9] = lookcount #Index in Zdata
        Zi[lookcount,10] = indexcount #Index Z
        Zi[lookcount,11] = i #Index in Data
        # println(size(data))
        for j in 1:nCh
            # println(typeind+j-1)
            # data[i,typeind+j-1]
            Zi[lookcount,11+j] = data[i,typeind+j-1] #Index in Data
        end
        if (i == size(data)[1])
            Zi = Zi[1:lookcount, :]
            push!(data2d,Zi)
        end
        # if data[i,zind] == 160
        #     println("Meto esto ", Zi[lookcount,5])
        # end
    end
    # indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5; indexA = 6; indexB = 7; indexTag = 8; indexData = 9; indexInZ = 10; indexInRawData = 11;
    data2d, data, numberObj2D
end

struct zParam
    minZ
    medZ
    maxZ
end
#cuando ni procesas y si al watersheed, vas a por todas.
# usadas por artificiales y figura 1 (reales de figura 1)
function ZmedReftoXY(folder::String; medZ = 0)
    data = CSV.read(folder, delim='\t', DataFrame)
    #
    #
    indexX = 3; indexY = 4; indexArea = 2; indexZ = 1; indexAng = 5; indexA = 8; indexB = 9;
    indexFeret = 10; indexFeretX = 11; indexFerety = 12; indexFeretAngle = 13; indexFeretmin = 14;
    indexBX = 15; indexBY= 16;indexArea = 2; indexBW= 15; indexBH= 16;
    abW = data[:,indexBW]; bbW  = data[:,indexBH];

    nSlices = data[end,indexZ]

    resliceFactor = 150/nSlices

    areaW =  (data[:, indexA]./2) .* (data[:, indexB]./2) .* pi;


    medArea = mean(areaW);
    var2 = [areaW[i] for i in 1:length(areaW) if medArea<areaW[i]];
    maxZ = sqrt(mean(var2)/pi)
    # var2 = [bbW[i] for i in 1:length(bbW) if value1<bbW[i]]; medZ = mean(var2);
    # medZ = (round(medZ));
    var2 = [areaW[i] for i in 1:length(areaW) if medArea>areaW[i]];
    minZ = sqrt(mean(var2)/pi)
    medZ = sqrt(medArea/pi)
    minZ = Int(round(minZ/resliceFactor))
    medZ = Int(round(medZ/resliceFactor))
    maxZ = Int.(round(maxZ/resliceFactor))
    if minZ < 2
        minZ = 2
    end
    return (zParam(minZ, medZ, maxZ))
end

function ZmedRefToPixelSize(folder::String, minR::Float64, pixWidth::Float64, voxDepth::Float64; medZ = 0)
    data = CSV.read(folder, delim='\t', DataFrame)
    #
    #
    indexX = 3; indexY = 4; indexArea = 2; indexZ = 1; indexAng = 5; indexA = 8; indexB = 9;
    indexFeret = 10; indexFeretX = 11; indexFerety = 12; indexFeretAngle = 13; indexFeretmin = 14;
    indexBX = 15; indexBY= 16;indexArea = 2; indexBW= 15; indexBH= 16;
    abW = data[:,indexBW]; bbW  = data[:,indexBH];


    areaW =  (data[:, indexA]./2) .* (data[:, indexB]./2) .* pi;

    medArea = minR^2*10
    medianPixelDiameter = sqrt(medArea/pi)*2
    medianMicraDiameterXY = medianPixelDiameter * pixWidth;
    medZ = Int(round(medianMicraDiameterXY/voxDepth));


    var2 = [areaW[i] for i in 1:length(areaW) if medArea<areaW[i]];
    maxArea = mean(var2)
    maxPixelDiameter = sqrt(maxArea/pi)*2
    maxMicraDiameterXY = maxPixelDiameter * pixWidth;
    maxZ = Int(round(maxMicraDiameterXY/voxDepth));

    var2 = [areaW[i] for i in 1:length(areaW) if medArea>areaW[i]];
    minArea = mean(var2)/pi
    minPixelDiameter = sqrt(minArea/pi)*2
    minMicraDiameterXY = minPixelDiameter * pixWidth;
    minZ = Int(round(minMicraDiameterXY/voxDepth));
    return (zParam(minZ, medZ, maxZ))
end

function Zmin_Zmed_Zmax(folder::String; medZ=0)
    #
    #
    data = CSV.read(folder, delim='\t')
    #
    #
    indexX = 3; indexY = 4; indexArea = 2; indexZ = 1; indexAng = 5; indexA = 6; indexB = 7;
    indexFeret = 8; indexFeretX = 9; indexFerety = 10; indexFeretAngle = 11; indexFeretmin = 12;
    indexBX = 13; indexBY= 14; indexArea = 2; indexBW= 15; indexBH= 16;
    # abW = data[:,indexBW];
    bbW_f  = data[:,indexBH];

    varW = data[:,indexArea]; f_AR_m = mean(varW); f_AR_s = std(varW);
    bbW = Vector{Float64}(); Ob_areaW  = Vector{Float64}();
    for i in 1:size(data)[1]
        if (f_AR_m - f_AR_s <= varW[i] >= f_AR_m - f_AR_s)
            push!(bbW, data[i,indexBH]);
            push!(Ob_areaW, data[i, indexArea]);
        end
    end
    # Ob_areaW = data[:,indexArea];
    # El_areaW =  (data[:, indexA]./2) .* (data[:, indexB]./2) .* pi;
    # shapeW = Ob_areaW./El_areaW; aspctRatW = data[:, indexB] ./ data[:, indexA];
    if medZ > 0
        var2 = [bbW[i] for i in 1:length(bbW) if medZ<bbW[i]];
        if length(var2) > 0
            maxZ = mean(var2);
        else
            maxZ = medZ;
        end
        var2 = [bbW[i] for i in 1:length(bbW) if medZ>bbW[i]];
        if length(var2) > 0
            minZ = mean(var2);
        else
            minZ = medZ;
        end
    else
        medZ = mean(bbW);
        #
        # minZ = medZ - Int.(round.(std(newV)));
        fThrs = medZ;
        var3 = [bbW[i] for i in 1:length(bbW) if (bbW[i] < medZ)];
        # var3 = [newV[i] for i in 1:length(newV) if (newV[i]<fThrs)];
        if length(var3) > 0
            minZ = mean(var3);
        else
            minZ = medZ;
        end
        #
        # maxZ = medZ + Int.(round.(std(newV)));
        fThrs = medZ;
        var3 = [bbW[i] for i in 1:length(bbW) if (bbW[i] > medZ)];
        # var3 = [newV[i] for i in 1:length(newV) if (newV[i]<fThrs)];
        if length(var3) > 0
            maxZ = mean(var3);
        else
            maxZ = medZ;
        end
        #
    end
    #
    minZ = Int(round(minZ));  medZ = Int(round(medZ));
    maxZ = Int.(round(maxZ));
    return (zParam(minZ, medZ, maxZ))
end
#esta es laque vale cuando procesas y watersheed en XZ
#usadas por el resto porque es un improve que creemos que va mejor para nuestra muestra (phenomoligcaly adjusted)
function min_max_Zdef(folder::String; medZ=0)
    #
    #
    data = CSV.read(folder, delim='\t', DataFrame)
    #
    #
    indexX = 3; indexY = 4; indexArea = 2; indexZ = 1; indexAng = 5; indexA = 6; indexB = 7;
    indexFeret = 8; indexFeretX = 9; indexFerety = 10; indexFeretAngle = 11; indexFeretmin = 12;
    indexBX = 13; indexBY= 14;indexArea = 2; indexBW= 15; indexBH= 16;
    abW = data[:,indexBW]; bbW  = data[:,indexBH];

    areaW =  (data[:, indexA]./2) .* (data[:, indexB]./2) .* pi;
    if medZ > 0
        var2 = [bbW[i] for i in 1:length(bbW) if medZ<bbW[i]];
        if length(var2) > 0
            maxZ = mean(var2);
        else
            maxZ = medZ;
        end
        var2 = [bbW[i] for i in 1:length(bbW) if medZ>bbW[i]];
        if length(var2) > 0
            minZ = mean(var2);
        else
            minZ = medZ;
        end
    else
        # indV = [i for i in 1:length(bbW) if percentile(bbW, 25)<bbW[i]<percentile(bbW,75)];
        # areaW = areaW[indV]; bbW = bbW[indV];
        # value1 = mean(areaW);
        # var2 = [bbW[i] for i in 1:length(areaW) if value1<areaW[i]]; medZ = mean(var2);
        value1 = mean(areaW);
        var2 = [bbW[i] for i in 1:length(areaW) if value1<areaW[i]];
        value1 = mean(var2);  medZ = (round(value1));
        # var2 = [bbW[i] for i in 1:length(bbW) if value1<bbW[i]]; medZ = mean(var2);
        # medZ = (round(medZ));
        var2 = [bbW[i] for i in 1:length(bbW) if medZ<bbW[i]];
        if isempty(var2) == false
            maxZ = mean(var2);
        else
            maxZ = medZ
        end
        # maxZ = mean(var2);
        var2 = [bbW[i] for i in 1:length(bbW) if medZ>bbW[i]];
        if isempty(var2) == false
            minZ = mean(var2);
        else
            minZ = medZ
        end
    end
    #
    minZ = Int(round(minZ));
    medZ = Int(round(medZ));
    maxZ = Int.(round(maxZ));
    return (zParam(minZ, medZ, maxZ))
end
function min_max_Zdef2(folder::String; medZ=0)
    #
    #
    data = CSV.read(folder, delim='\t', DataFrame)
    #
    #
    indexX = 3; indexY = 4; indexArea = 2; indexZ = 1; indexAng = 5; indexA = 6; indexB = 7;
    indexFeret = 8; indexFeretX = 9; indexFerety = 10; indexFeretAngle = 11; indexFeretmin = 12;
    indexBX = 13; indexBY= 14;indexArea = 2; indexBW= 15; indexBH= 16;
    abW = data[:,indexBW]; bbW  = data[:,indexBH];

    areaW =  (data[:, indexA]./2) .* (data[:, indexB]./2) .* pi;
    if medZ > 0
        var2 = [bbW[i] for i in 1:length(bbW) if medZ<bbW[i]];
        if length(var2) > 0
            maxZ = mean(var2);
        else
            maxZ = medZ;
        end
        var2 = [bbW[i] for i in 1:length(bbW) if medZ>bbW[i]];
        if length(var2) > 0
            minZ = mean(var2);
        else
            minZ = medZ;
        end
    else
        # indV = [i for i in 1:length(bbW) if percentile(bbW, 25)<bbW[i]<percentile(bbW,75)];
        # areaW = areaW[indV]; bbW = bbW[indV];
        # value1 = mean(areaW);
        # var2 = [bbW[i] for i in 1:length(areaW) if value1<areaW[i]]; medZ = mean(var2);
        value1 = mean(areaW);
        var2 = [bbW[i] for i in 1:length(areaW) if value1<areaW[i]];
        value1 = mean(var2);  medZ = (round(bbW));
        # var2 = [bbW[i] for i in 1:length(bbW) if value1<bbW[i]]; medZ = mean(var2);
        # medZ = (round(medZ));
        var2 = [bbW[i] for i in 1:length(bbW) if medZ<bbW[i]];
        if isempty(var2) == false
            maxZ = mean(var2);
        else
            maxZ = medZ
        end
        # maxZ = mean(var2);
        var2 = [bbW[i] for i in 1:length(bbW) if medZ>bbW[i]];
        if isempty(var2) == false
            minZ = mean(var2);
        else
            minZ = medZ
        end
    end
    #
    minZ = Int(round(minZ));
    medZ = Int(round(medZ));
    maxZ = Int.(round(maxZ));
    return (zParam(minZ, medZ, maxZ))
end
#esta va a chorro, priemro min y luego med. Va mejor con XZ procesado (juntito)
function Zmin_Zmed_Zmax2(folder::String, medZ=0)
    #
    #
    data = CSV.read(folder, delim='\t')
    #
    #
    indexX = 3; indexY = 4; indexArea = 2; indexZ = 1; indexAng = 5; indexA = 6; indexB = 7;
    indexFeret = 8; indexFeretX = 9; indexFerety = 10; indexFeretAngle = 11; indexFeretmin = 12;
    indexBX = 13; indexBY= 14; indexArea = 2; indexBW= 15; indexBH= 16;
    # abW = data[:,indexBW];
    bbW_f  = data[:,indexBH];
    if medZ > 0
        var2 = [bbW[i] for i in 1:length(bbW) if medZ<bbW[i]];
        if length(var2) > 0
            maxZ = mean(var2);
        else
            maxZ = medZ;
        end
        var2 = [bbW[i] for i in 1:length(bbW) if medZ>bbW[i]];
        if length(var2) > 0
            minZ = mean(var2);
        else
            minZ = medZ;
        end
    else
        varW = data[:,indexArea]; f_AR_m = mean(varW); f_AR_s = std(varW);
        bbW1 = Vector{Float64}(); Ob_areaW1  = Vector{Float64}();
        bbW2 = Vector{Float64}(); Ob_areaW2  = Vector{Float64}();
        for i in 1:size(data)[1]
            if (varW[i] < f_AR_m)
                push!(bbW1, data[i,indexBH]);
                push!(Ob_areaW1, data[i, indexArea]);
            else
                push!(bbW2, data[i,indexBH]);
                push!(Ob_areaW2, data[i, indexArea]);
            end
        end
        minZ = (mean(bbW1)); maxZ = (mean(bbW2));
        #
        var3 = [bbW_f[i] for i in 1:length(bbW_f) if (maxZ > bbW_f[i] > minZ)];
        medZ = mean(var3);
        #
    end
    #
    minZ = Int(round(minZ));  medZ = Int(round(medZ));
    maxZ = Int.(round(maxZ));
    return (zParam(minZ, medZ, maxZ))
end

function sortPoints_ascending(xref,yref,xar,yar)
    dist = zeros(length(xar))
    for i in 1:length(xar)
         dist[i] = sqrt((xar[i].-xref).^2 .+ (yar[i].-yref).^2)
        #dist[i] = ((xar[i].-xref).^2 .+ (yar[i].-yref).^2)
    end
    inMinSorted = sortperm(vec(dist))
    B = dist[inMinSorted]
    #dist = B
    return inMinSorted, B, dist
end

function generate_points(xref,yref,angref,aref,bref)
    # This function works with a clockwise rotation matrix, so values of angles
    # from IJ (where the Y axis is inverted) needs to be multiplied by -1 when
    # they are used as function arguments. Moreover, angles values must be
    # defined in angle degrees
    angref = (angref*pi)/180; nn = 359;
    theta = range(0, 2*pi, length=nn); coordEl = zeros(2, nn);
    coordEl[1, :] = aref .* cos.(theta); coordEl[2, :] = bref .* sin.(theta);

    rot_matrix = [cos(angref) -sin(angref); sin(angref) cos(angref)];
    coordN = rot_matrix * coordEl;
    xpins = coordN[1, :] .+ xref; ypins = coordN[2, :] .+ yref;
    return(xpins, ypins);
end

struct overlappedInfo
    xcoordA
    ycoordA
    xcoordB
    ycoordB
    AonB
    BonA
end
function OverlappedEllipses(aA, bA, xA, yA, angA, aB, bB, xB, yB, angB)
     # A on B
    xdA, ydA = generate_points(xA, yA, angA, aA, bA)
    angB = (angB*pi)/180; cB = sqrt(abs((aB)^2-(bB)^2));
    f_1_x_B = xB-cB*cos(angB); f_1_y_B = yB-cB*sin(angB);
    f_2_x_B = xB+cB*cos(angB); f_2_y_B = yB+cB*sin(angB);
    qA = 0
    for k in 1:length(xdA)
        dis1tB = sqrt((xdA[k].-f_1_x_B).^2 .+(ydA[k].-f_1_y_B).^2)
        dis2tB = sqrt((xdA[k].-f_2_x_B).^2 .+(ydA[k].-f_2_y_B).^2)
        if (dis1tB+dis2tB) <= 2*aB
              qA = qA + 1;
        end
    end
    # B on A
    xdB, ydB = generate_points(xB, yB, angB, aB, bB)
    angA = (angA*pi)/180; cA = sqrt(abs((aA)^2-(bA)^2));
    f_1_x_A = xA-cA*cos(angA); f_1_y_A = yA-cA*sin(angA);
    f_2_x_A = xA+cA*cos(angA); f_2_y_A = yA+cA*sin(angA);
    qB = 0
    for k in 1:length(xdB)
        dis1tA = sqrt((xdB[k].-f_1_x_A).^2 .+(ydB[k].-f_1_y_A).^2)
        dis2tA = sqrt((xdB[k].-f_2_x_A).^2 .+(ydB[k].-f_2_y_A).^2)
        if (dis1tA+dis2tA) <= 2*aA
              qB = qB + 1;
        end
    end
    overlappedInfo(xdA, ydA, xdB, ydB, qA, qB)
end

struct Line
    beta
    ypredicted
    y
    supCI
    downCI
    Rsquared
end
# 1 arg normal, 1 arg opcional, 1 argumento que no es opcional, pero que no imoprta donde lo pongas
function LinearRegression(y, x::Array{Float64} = hcat(ones(length(y)),1:1:length(y)); points_to_fit::Int64 = length(y))
    if points_to_fit != length(y)
        y1 = y[1:points_to_fit]
        y2 = y
        x1 = x[1:points_to_fit,:]
        x2 = x
    else
        y1 = y2 = y
        x1 = x2 = x
    end
    beta = inv(x1'*x1)*(x1'*y1)
    error = y1 - x1*beta
    sigmaSquared = Statistics.var(error)*(length(y1)-1)/(length(y1)-size(beta,2))
    sigma = sqrt(sigmaSquared)*1.86
    y2 = y
    ypredicted = x2*beta
    supCI = ypredicted .+ sigma
    downCI = ypredicted .- sigma
    SSerr = sum((y2 .- ypredicted).^2)
    SStot = sum((y2 .- mean(y2)).^2)
    #MSE = SE/length(y2)
    # RMSE = sqrt(MSE)
    #rMSE = MSE/var(y2)
    #Rsquared = 1 - rMSE
    Rsquared = 1 - SSerr/SStot
    return Line(beta, ypredicted, y2, supCI, downCI, Rsquared)
end

struct Vector3D
    xpred
    ypred
    xcenter
    ycenter
    zcenter
    xref    # punto prdicho en el plano siguiente
    yref
    betaX
    betaY
    vx
    vy
    vz
    vmodule
end
# todos los argumentos son opcionales, y van en orden
function Line_3D(;x::Array{Float64}=ones(0), y::Array{Float64}=ones(0), z::Array{Float64}=collect(range(1.0,length(x),step =1)), points_to_fit::Int64 = length(x))
    #collect(range(1, stop=length(x), length=length(x)))
    #println(z);
    z = hcat(ones(length(z)),z);
    regX = LinearRegression(x,z, points_to_fit= points_to_fit)
    regY = LinearRegression(y,z, points_to_fit= points_to_fit)
    xpred = regX.ypredicted
    ypred = regY.ypredicted
    xcenter = (xpred[1]+xpred[end])/2
    ycenter = (ypred[1]+ypred[end])/2
    zcenter = (z[1,2]+z[end,2])/2
    betaX = regX.beta; betaY = regY.beta;
    # betaX = LinearRegression(x,z, points_to_fit=points_to_fit).beta #No haria falta hacerlo otra vez no?
    # betaY = LinearRegression(y,z, points_to_fit=points_to_fit).beta
    # xref = xpred[end] #[1 (z[end,2]+1)] * betaX
    # yref = ypred[end] #[1 (z[end,2]+1)] * betaY
    # # es lo mismo en matrizx que en lo otro, `pero el output te da como array y el Float64() solo acepta int
    # xref = xref[1]; yref = yref[1];
    xref = betaX[2] * z[end,2]+1 + betaX[1];
    yref = betaY[2] * z[end,2]+1 + betaY[1];
    vx = xpred[end] - xpred[1]
    vy = ypred[end] - ypred[1]
    vz = z[end,2] - z[1,2]
    vmodule = sqrt(vx^2+vy^2+vz^2)
    return Vector3D(xpred, ypred, xcenter, ycenter, zcenter, xref, yref, betaX, betaY, vx, vy, vz, vmodule)
end

function outliersDetection(vDistV, sliceVector, medSizeMin)
    #ztest
    a = vDistV[1:medSizeMin];
    ztest = hcat(ones(length(vDistV[1:medSizeMin])), range(1, length(vDistV[1:medSizeMin]), step=1))
    if det(ztest'*ztest) != 0
        regD = LinearRegression(vDistV, points_to_fit = medSizeMin)
        ypredicted = regD.ypredicted; y = regD.y; supCI = regD.supCI;
        # downCI = regD.downCI; Rsquared = regD.Rsquared;
        # sigma = findmax(y[1:medSizeMin])[2][1] - findmin(y[1:medSizeMin])[2][1]
        sigma = abs(percentile(y[1:medSizeMin], 25) - percentile(y[1:medSizeMin], 75));
        # sigma = percentile(y[1:medSizeMin], 0) - percentile(y[1:medSizeMin], 100);
        supCI = ypredicted .+ sigma
        downCI = ypredicted .- sigma
    else
        y = vDistV;
        # sigma = findmax(y[1:medSizeMin])[2][1] - findmin(y[1:medSizeMin])[2][1]
        sigma = percentile(y[1:medSizeMin], 25) - percentile(y[1:medSizeMin], 75);
        ypredicted = ones(length(y), 1) .* median(y[1:medSizeMin]);
        supCI = ypredicted .+ sigma
        downCI = ypredicted .- sigma
    end
    #
    H = sliceVector; y1 = Vector{Float64}();
    for point in 1:length(y)
        #b = findmin(abs.(ypredicted.-y[point]))[2]
        b = point
        if y[point] > supCI[b] || y[point] < downCI[b]
            H[point] = H[point].*-1
        else
            push!(y1,vDistV[point])
        end
    end
    cutPoints = Vector{Int64}(); returnPoints = Vector{Int64}();
    for i in 1:length(H)
        if H[i] < 0
            push!(cutPoints, i);  c = 1;
            for j in i:length(H)-1
                if H[j+1] > 0
                    c = 0; push!(returnPoints, j)
                    break;
                end
            end
            if c == 1
                push!(returnPoints, length(H))
            end
        end
    end
    H, y1, cutPoints, returnPoints
end

function volDelimiter(sizeVec, minVol, maxVol)
    cutPoint_min = 0; cutPoint_max = length(sizeVec);
    for imaxt in length(sizeVec):-1:1
        if sizeVec[imaxt] >= maxVol
            cutPoint_max = imaxt;
            maxVol = sizeVec[end];
        elseif (sizeVec[imaxt]) <= minVol
            cutPoint_min = imaxt;
            break;
        end
    end
    return(cutPoint_min, cutPoint_max)
end

function angPreobj(xx, yy, zz, testV_pre)
    angV = Vector{Float64}();
    # vector reference
    ztest = hcat(ones(length(xx)), zz)
    if det(ztest'*ztest) != 0
        line = Line_3D(x = xx, y = yy, z = zz)
        a = line.betaX; pX = a[2] #slope
        b = line.betaY; pY = b[2] #slope
    else
        pX = 0; pY = 0;
    end
    refV = [pX pY (zz[end] - zz[1])]
    #normalizar ambos vectores para no cagarla
    # Angles between the reference vector and the list of vectors
    for i in 1:size(testV_pre, 1)
        a = testV_pre[i, 1]; b = testV_pre[i, 2];
        testV = [a b 1]
        angt = cos(dot(testV, refV) / (norm(testV) * norm(refV)))
        push!(angV, angt)
    end
    return(angV, refV)
end

function getCovariance(var1, var2)
    pearson_ind = nothing;
    if (std(var1) * std(var2)) == 0
        pearson_ind = 1;
    else
        pearson_ind = cov(var1, var2) / (std(var1) * std(var2))
    end
    pearson_ind
end

struct infoPreobj
    dist
    ang
    pearson
end
function info_Preobj(tempEl, index1=1, index2=size(tempEl, 1))
    indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5; indexA = 6; indexB = 7; indexTag = 8; indexData = 9; indexInZ = 10; indexInRawData = 11;
    #
    if index1 < 1 || index2 > size(tempEl, 1)
        xS = tempEl[:, indexX]
        yS = tempEl[:, indexY]
        zS = tempEl[:, indexZ]
        n_points = Int(length(xS));
    else
        xS = tempEl[index1:index2, indexX]
        yS = tempEl[index1:index2, indexY]
        zS = tempEl[index1:index2, indexZ]
        n_points = Int(length(xS));
    end
    ztest = hcat(ones(length(xS)), collect(range(1.0,length(xS),step =1)));
    if det(ztest'*ztest) != 0
        # line = Line_3D(x = xS, y = yS, points_to_fit = n_points)
        # println("zS = ", zS)
        line = Line_3D(x = xS, y = yS, z = zS)
        cx = line.xcenter; cy = line.ycenter;
        a = line.betaX; pX = a[2] #slope
        b = line.betaY; pY = b[2] #slope
    else
        cx = median(xS); cy = median(yS);
        pX = 0; pY = 0;
    end
    distV = Vector{Float64}(); angV = Vector{Float64}(); pearson_ind = 0;
    if size(tempEl, 1) > 1
        #distance
        indV_sorted, distV_sorted, distV = sortPoints_ascending(cx, cy, tempEl[:,indexX], tempEl[:,indexY])
        #angle
        refV = [pX pY (tempEl[end, indexZ] - tempEl[1, indexZ])]
        for i in 1:size(tempEl, 1)-1
            dx = tempEl[i+1, indexX] - tempEl[i, indexX]
            dy = tempEl[i+1, indexY] - tempEl[i, indexY]
            testV = [dx dy 1]
            angt = cos(dot(testV, refV) / (norm(testV) * norm(refV)))
            push!(angV, angt)
        end
        #pearson
        # pearsonV = Vector{Float64}();
        if length(distV) > 1
            a = distV[2:end]; b = angV;
            pearson_ind = getCovariance(distV[2:end], angV);
        end
        # push!(pearsonV, pearson_ind)
    end
    infoPreobj(distV, angV, pearson_ind)
end

function info_Preobj_distZ(tempEl, index1=1, index2=size(tempEl, 1))
    indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5; indexA = 6; indexB = 7; indexTag = 8; indexData = 9; indexInZ = 10; indexInRawData = 11;
    #
    if index1 < 1 || index2 > size(tempEl, 1)
        xS = tempEl[:, indexX]
        yS = tempEl[:, indexY]
        zS = tempEl[:, indexZ]
        n_points = Int(length(xS));
    else
        xS = tempEl[index1:index2, indexX]
        yS = tempEl[index1:index2, indexY]
        zS = tempEl[index1:index2, indexZ]
        n_points = Int(length(xS));
    end
    ztest = hcat(ones(length(xS)), collect(range(1.0,length(xS),step =1)));
    if det(ztest'*ztest) != 0
        # line = Line_3D(x = xS, y = yS, points_to_fit = n_points)
        # println("zS = ", zS)
        line = Line_3D(x = xS, y = yS, z = zS)
        cx = line.xcenter; cy = line.ycenter;
        a = line.betaX; pX = a[2] #slope
        b = line.betaY; pY = b[2] #slope
        xpred = hcat(ones(length(tempEl[:,indexZ])),tempEl[:,indexZ])*a;
        ypred = hcat(ones(length(tempEl[:,indexZ])),tempEl[:,indexZ])*b;
    else
        cx = median(xS); cy = median(yS);
        pX = 0; pY = 0;
    end
    distV = Vector{Float64}(); angV = Vector{Float64}(); pearson_ind = 0;
    if size(tempEl, 1) > 1
        #distance
        try
            distV = sqrt.((xpred.-tempEl[:,indexX]).^2 .+ (ypred .- tempEl[:,indexY]) .^2)
        catch e
            println("No predicted line")
            indV_sorted, distV_sorted, distV = sortPoints_ascending(cx, cy, tempEl[:,indexX], tempEl[:,indexY])
        end
        #angle
        refV = [pX pY (tempEl[end, indexZ] - tempEl[1, indexZ])]
        for i in 1:size(tempEl, 1)-1
            dx = tempEl[i+1, indexX] - tempEl[i, indexX]
            dy = tempEl[i+1, indexY] - tempEl[i, indexY]
            testV = [dx dy 1]
            angt = cos(dot(testV, refV) / (norm(testV) * norm(refV)))
            push!(angV, angt)
        end
        #pearson
        # pearsonV = Vector{Float64}();
        if length(distV) > 1
            a = distV[2:end]; b = angV;
            pearson_ind = getCovariance(distV[2:end], angV);
        end
        # push!(pearsonV, pearson_ind)
    end
    infoPreobj(distV, angV, pearson_ind)
end

function summary_objs(RefStruct, nCh)
    indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5; indexA = 6; indexB = 7; indexTag = 8; indexData = 9; indexInZ = 10; indexInRawData = 11;
    indexType = 12;
    objects = RefStruct.preok;
    # objects3D_df = DataFrame(x = Float64[], y = Float64[], z = Float64[], angXY = Float64[], A = Float64[], B = Float64[], C = Float64[], vol = Float64[], vx = Float64[], vy = Float64[], vz = Float64[])
    matrixData = zeros(size(objects,1), indexType+nCh-1); c = 0;
    for i in 1:size(objects,1)
        x = objects[i][:,indexX]
        y = objects[i][:,indexY]
        z = objects[i][:,indexZ]
        angle = objects[i][:,indexAng]
        a = objects[i][:,indexA]
        b = objects[i][:,indexB]
        area = objects[i][:,indexArea]
        ztest = hcat(ones(length(x)), collect(range(1.0,length(x),step =1)));
        if det(ztest'*ztest) != 0
            # line = Line_3D(x = xS, y = yS, points_to_fit = n_points)
            vector3D = Line_3D(x = x, y = y, z = z)
            cx = vector3D.xcenter; cy = vector3D.ycenter;
            cz = vector3D.zcenter;
            vx = vector3D.vx; vy = vector3D.vy; vz = vector3D.vz;
        else
            cx = median(x); cy = median(y);
            cz = median(z);
            vx = 0.0; vy = 0.0; vz = 0.0;
        end
        medianA = median(a)
        medianB = median(b)
        medianC = medianC = sqrt(vx^2+vy^2+vz^2)
        #z[end] - z[1];
        medianAngleXY = median(angle)
        vol = sum(area);
        # push!(objects3D_df, (cx, cy, cz, medianAngleXY, medianA, medianB, medianC, vol, vx, vy, vz))
        #
        dataV = hcat(cx, cy, cz, medianAngleXY, medianA, medianB, medianC, vol, vx, vy, vz)
        for j in 1:nCh
            tempt = (objects[i][:,indexType+j-1])
            dataV = hcat(dataV, sum(tempt));
        end
        # for j in 1:nCh
        #     tempt = Int.(objects[i][:,indexType+j-1])
        #     dataV = hcat(dataV, mean(tempt));
        # end
        #
        # push!()
        c = c + 1; matrixData[c, :] = dataV;
    end
    # objects3D_df = convert(DataFrame, matrixData)
    names = vcat([:Xcenter, :Ycenter, :Zcenter, :AngleXY, :A, :B, :C, :Vol, :Vx, :Vy, :Vz], Symbol.(:Channel,axes(ones(nCh),1)))
    objects3D_df = DataFrame(Tables.table(matrixData, header = names))

    # for j in 1:nCh
    #     names = vcat(names, Symbol("Channel$j"))
    # end
    # objects3D_df = DataFrames.rename(objects3D_df, names)

    # objects3D = DataFrame(x = vector3D.xcenter, y = vector3D.ycenter, z = vector3D.zcenter, angXY = medianAngleXY, A = medianA, B = medianB, C = medianC, TypeMed = medianType, TypeMode = modeType, vx = vector3D.vx, vy = vector3D.vy, vz =  vector3D.vz)
    # objects3D = vcat(objects3D, hcat(vector3D.xcenter,vector3D.ycenter, vector3D.zcenter, vector3D.vx, vector3D.vy, vector3D.vz, vector3D.vmodule, medianType, modeType))
    return(objects3D_df)
end
function summary_from_folder(folder, nCh)
    dir = readdir(folder)
    preok = Vector{}()
    for i in 1:length(dir)
        tF = dir[i];
        if tF[1:6] == "Object"
            println(tF);
            push!(preok, readdlm(folder*"\\"*tF, '\t', header = true)[1])
        end
    end
    indX = 1; indY = 2; indZ = 3; indAr = 4; indAng = 5; indA = 6; indB = 7;
    data2d_index = 8; Zindex_data2d = 9; rawData_index = 10;
    indexCh_in = 11;
    objects = preok;
    # objects3D_df = DataFrame(x = Float64[], y = Float64[], z = Float64[], angXY = Float64[], A = Float64[], B = Float64[], C = Float64[], vol = Float64[], vx = Float64[], vy = Float64[], vz = Float64[])
    matrixData = zeros(size(objects,1), indexCh_in+nCh); c = 0;
    for i in 1:size(objects,1)
        x = objects[i][:,indX]
        y = objects[i][:,indY]
        z = objects[i][:,indZ]
        angle = objects[i][:,indAng]
        a = objects[i][:,indA]
        b = objects[i][:,indB]
        area = objects[i][:,indAr]
        ztest = hcat(ones(length(x)), collect(range(1.0,length(x),step =1)));
        if det(ztest'*ztest) != 0
            # line = Line_3D(x = xS, y = yS, points_to_fit = n_points)
            vector3D = Line_3D(x = x, y = y, z = z)
            cx = vector3D.xcenter; cy = vector3D.ycenter;
            cz = vector3D.zcenter;
            vx = vector3D.vx; vy = vector3D.vy; vz = vector3D.vz;
        else
            cx = median(x); cy = median(y);
            cz = median(z);
            vx = 0.0; vy = 0.0; vz = 0.0;
        end
        medianA = median(a)
        medianB = median(b)
        medianC = medianC = sqrt(vx^2+vy^2+vz^2)
        #z[end] - z[1];
        medianAngleXY = median(angle)
        vol = sum(area);
        # push!(objects3D_df, (cx, cy, cz, medianAngleXY, medianA, medianB, medianC, vol, vx, vy, vz))
        #
        dataV = hcat(cx, cy, cz, medianAngleXY, medianA, medianB, medianC, vol, vx, vy, vz)
        for j in 1:nCh
            tempt = Int.(objects[i][:,indexCh_in+j-1])
            dataV = hcat(dataV, sum(tempt));
        end
        # for j in 1:nCh
        #     tempt = Int.(objects[i][:,indexCh_in+j-1])
        #     dataV = hcat(dataV, mean(tempt));
        # end
        #
        # push!()
        c = c + 1; matrixData[c, :] = dataV;
    end
    objects3D_df = convert(DataFrame, matrixData)
    names = [:Xcenter, :Ycenter, :Zcenter, :AngleXY, :A, :B, :C, :Vol, :Vx, :Vy, :Vz]
    for j in 1:nCh
        names = vcat(names, Symbol("Channel$j"))
    end
    names!(objects3D_df, names)
    # objects3D = DataFrame(x = vector3D.xcenter, y = vector3D.ycenter, z = vector3D.zcenter, angXY = medianAngleXY, A = medianA, B = medianB, C = medianC, TypeMed = medianType, TypeMode = modeType, vx = vector3D.vx, vy = vector3D.vy, vz =  vector3D.vz)
    # objects3D = vcat(objects3D, hcat(vector3D.xcenter,vector3D.ycenter, vector3D.zcenter, vector3D.vx, vector3D.vy, vector3D.vz, vector3D.vmodule, medianType, modeType))
    return(objects3D_df)
end
function objectToRow_3D(obj3d, nCh)
    #index for object
    indX = 1; indY = 2; indZ = 3; indAr = 4; indAng = 5; indA = 6; indB = 7;
    data2d_index = 8; Zindex_data2d = 9; rawData_index = 10;
    indexCh_in = 11;
    #index for Line: this is the column indexing for each variable
    # indexType = 12;
    # indexX = 1; indexY = 2; indexZ = 3; indexArea = 4; indexAng = 5; indexA = 6; indexB = 7;
    # indxVol = 8; indexVx = 9; indexVy = 10; indexVz = 11;
    x = obj3d[:,indX]
    y = obj3d[:,indY]
    z = obj3d[:,indZ]
    angle = obj3d[:,indAng]
    a = obj3d[:,indA]
    b = obj3d[:,indB]
    area = obj3d[:,indAr]
    ztest = hcat(ones(length(x)), collect(range(1.0,length(x),step =1)));
    if det(ztest'*ztest) != 0
        # line = Line_3D(x = xS, y = yS, points_to_fit = n_points)
        vector3D = Line_3D(x = x, y = y, z = z)
        cx = vector3D.xcenter; cy = vector3D.ycenter;
        cz = vector3D.zcenter;
        vx = vector3D.vx; vy = vector3D.vy; vz = vector3D.vz;
    else
        cx = median(x); cy = median(y);
        cz = median(z);
        vx = 0.0; vy = 0.0; vz = 0.0;
    end
    medianA = median(a)
    medianB = median(b)
    medianC = sqrt(vx^2+vy^2+vz^2)
    #z[end] - z[1];
    medianAngleXY = median(angle)
    vol = sum(area);
    # push!(objects3D_df, (cx, cy, cz, medianAngleXY, medianA, medianB, medianC, vol, vx, vy, vz))
    #
    dataV = hcat(cx, cy, cz, medianAngleXY, medianA, medianB, medianC, vol, vx, vy, vz)
    # println(size(obj3d));
    for j in 1:nCh
        # println(indexCh_in+j-1)
        tempt = (obj3d[:,indexCh_in+j-1])
        dataV = hcat(dataV, sum(tempt));
    end
    dataV
end

function toDataFrame(RefStruct, dir, name, nCh)
   indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5; indexA = 6; indexB = 7; indexTag = 8; indexData2d = 9; indexInZ = 10; indexInRawData = 11;
   indexType = 12;
   objects = RefStruct.preok;
   obj = nothing;
   objects3D = zeros(1,9);
   for i in 1:size(objects,1)
       x = objects[i][:,indexX]; y = objects[i][:,indexY]; z = objects[i][:,indexZ];
       area = objects[i][:,indexA].*objects[i][:,indexB].*pi; angle = objects[i][:,indexAng]
       a = objects[i][:,indexA]; b = objects[i][:,indexB];
       indexdata2d = objects[i][:,indexData2d];
       indexZdata2d = objects[i][:,indexInZ];
       indexData = objects[i][:,indexInRawData];
       obj = hcat(x,y,z,area,angle,a,b,indexdata2d,indexData)
       for j in 1:nCh
           tempt = objects[i][:,indexType+j-1];
           obj = hcat(obj, tempt)
       end
       names = vcat([:Xcenter, :Ycenter, :Zcenter, :Area, :Angle, :A, :B, :IndexData2D, :IndexData], Symbol.(:Channel,axes(ones(nCh),1)))
       obj = DataFrame(Tables.table(obj, header = names))
       # obj = DataFrame(obj,:auto);
       # names!(obj,[:xc, :yc, :zc, :area, :angle, :majorc, :minorc, :type, :indData2D, :indZdata2D, :indexData])
       CSV.write("$dir\\Object$i"*"_$name", obj, delim = '\t')
       #name va ya con .txt
   end
   obj
end

struct jointEllipses_output
    nextEl_2
    nextStep
    cx # projected point to search nearest points by Line_3D
    cy
end
function EllipsesConnector(data2d_copy, z, tempEl, nextEl_1, medSizeMin)
    indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5; indexA = 6; indexB = 7; indexTag = 8; indexData = 9; indexInZ = 10; indexInRawData = 11;
    nextStep = 0; nextEl_2 = nothing;
    n_points = Int(round(medSizeMin));
    # projected cx and cy
    if size(tempEl,1) >= n_points
        # xS = tempEl[end-n_points+1:end, indexX]
        # yS = tempEl[end-n_points+1:end, indexY]
        # zS = tempEl[end-n_points+1:end, indexZ]
        # angref = median(tempEl[end-n_points+1:end, indexAng]);
        # aref = median(tempEl[end-n_points+1:end, indexA])
        # bref = median(tempEl[end-n_points+1:end, indexB]);
        xS = tempEl[1:n_points, indexX]
        yS = tempEl[1:n_points, indexY]
        zS = tempEl[1:n_points, indexZ]
        angref = median(tempEl[1:n_points, indexAng]);
        aref = median(tempEl[1:n_points, indexA])
        bref = median(tempEl[1:n_points, indexB]);
    else
        xS = tempEl[:, indexX]
        yS = tempEl[:, indexY]
        zS = tempEl[:, indexZ]
        angref = median(tempEl[:, indexAng]);
        aref = median(tempEl[:, indexA])
        bref = median(tempEl[:, indexB]);
    end
    ztest = hcat(ones(length(xS)), collect(range(1.0,length(xS),step =1)));
    if det(ztest'*ztest) != 0
        # line = Line_3D(x = xS, y = yS, points_to_fit = n_points)
        line = Line_3D(x = xS, y = yS, z = zS)
        cx = line.xref; cy = line.yref;
    else
        cx = median(xS); cy = median(yS);
    end
    #
    aup = tempEl[end, indexA]./2;
    bup = tempEl[end, indexB]./2;
    angup = tempEl[end, indexAng].*(-1);
    xup = tempEl[end, indexX];
    yup = tempEl[end, indexY]
    #
    indV, distV, dists = sortPoints_ascending(xup, yup, data2d_copy[z+1][:,indexX], data2d_copy[z+1][:,indexY])
    #
    k = 4
    if length(indV) < k
        k = length(indV)
    end
    in_sel = Vector{}();
    for i in 1:k
        at = (data2d_copy[z+1][indV[i],indexA]./2);
        bt = (data2d_copy[z+1][indV[i],indexB]./2);
        angt = data2d_copy[z+1][indV[i],indexAng].*(-1);
        xt = (data2d_copy[z+1][indV[i],indexX]);
        yt = (data2d_copy[z+1][indV[i],indexY]);
        #
        shared = OverlappedEllipses(aup, bup, xup, yup, angup, at, bt, xt, yt, angt)
        percOv_pair = shared.AonB + shared.BonA;
        #
        if percOv_pair > 0
            push!(in_sel, indV[i])
        else
            break
        end
    end
    #
    if length(in_sel) > 0
        indV2, distV2, dists2 = sortPoints_ascending(cx, cy, data2d_copy[z+1][in_sel, indexX], data2d_copy[z+1][in_sel,indexY]);
        #
        nextStep = 1;
        nextEl_2 = in_sel[indV2[1]];
        # at = (data2d_copy[z+1][nextEl_2,indexA]./2);
        # bt = (data2d_copy[z+1][nextEl_2,indexB]./2);
        # angt = data2d_copy[z+1][nextEl_2,indexAng].*(-1);
        # xt = (data2d_copy[z+1][nextEl_2,indexX]);
        # yt = (data2d_copy[z+1][nextEl_2,indexY]);
        # sharedL = OverlappedEllipses(aup, bup, xup, yup, angup, at, bt, xt, yt, angt)
        # percOv_pairL = sharedL.AonB + sharedL.BonA;
        # #
        # if percOv_pairL > 0
        #     nextStep = 1;
        # end
    end
    #
    jointEllipses_output(nextEl_2, nextStep, cx, cy)
end

struct preObj_output
    tempEl
    vVolV
    vDist
    vAng
    vPearson
end
function preObjElongation(sliceI, indexI, data2d, maxZinterval, minZinterval, medZinterval)
    data2d_copy = deepcopy(data2d)
    indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5; indexA = 6; indexB = 7; indexType = 8; indexTag = 8; indexData = 9; indexInZ = 10; indexInRawData = 11;
    # println("a ", data2d_copy[sliceI][indexI,12])
    # println("b ", data2d_copy[sliceI][indexI,13])
    # outputs def.
    vVolV = Vector{Float64}();
    # vDist = Vector{Float64}(); vAng = Vector{Float64}()
    # vPearson = nothing;
    # Initiator
    nextEl_1 = indexI; z = sliceI;
    tempEl = data2d_copy[z][nextEl_1,:]';
    pearson_ind = 0; volPre = tempEl[1, indexA] .* tempEl[1, indexB] .* pi;
    push!(vVolV, volPre);
    #elongation
    while z < sliceI+maxZinterval-1 && z < size(data2d_copy,1)
        ellipseInfo = EllipsesConnector(data2d_copy, z, tempEl, nextEl_1, minZinterval);
        goQ = ellipseInfo.nextStep; nextEl_2 = 0;
        if goQ == 1
            nextEl_2 = ellipseInfo.nextEl_2;
            #
            # println("nextEl_2 = ", nextEl_2)
            # println("previous Z = ", z)
            nextEl_1 = nextEl_2; z = z + 1;
            # println("nextEl_1 = ", nextEl_1)    #ojito que ya no es z+1
            # println("next Z = ", z)
            # println("que cojones pone en la línea? = ", data2d_copy[z][nextEl_1,:]')
            tempEl = vcat(tempEl, data2d_copy[z][nextEl_1,:]')
            #
            volPre = sum(tempEl[:, indexA] .* tempEl[:, indexB] .* pi)
            push!(vVolV, volPre);
            #
        else
            break
        end
    end
    #
    info = info_Preobj_distZ(tempEl, 1, medZinterval)
    distV = info.dist; angV = info.ang; pearson = info.pearson;
    vDist = distV; vAng = angV; vPearson = pearson;
    # println("a' ", tempEl[:,12])
    # println("b' ", tempEl[:,13])
    #
    preObj_output(tempEl, vVolV, vDist, vAng, vPearson)
end

struct TerminatorToOSCAR
    Hdist
    Hvect
    Htotal
    numbIn
    numbOut
    cutOut
    lengthOut
end
function TerminatorReturns(vDistV, angPreok, n_points)
    # println("size of d ", size(vDistV))
    # println("size of v ", size(angPreok))
    # println("n_points ", n_points)
    ##PROBANDO LAS FUNCIONES DE ARRIBA SON LAS BUENA outliersDetection
    HD, yD, cutPointsD, returnPointsD = outliersDetection(vDistV, collect(1:1:length(vDistV)), n_points)
    HV, yV, cutPointsV, returnPointsV = outliersDetection(angPreok, collect(1:1:length(angPreok)), n_points)
    # HD, yD, cutPointsD, returnPointsD = outliers_from_zero(info.dist[2:end])
    # HV, yV, cutPointsV, returnPointsV = outliers_from_zero(info.ang)
    # println("size of Hd ", size(HD))
    # println("size of Hv ", size(HV))
    HD = Int.(((HD ./ abs.(HD)) .- 1)./(-2));
    HV = Int.(((HV ./ abs.(HV)) .- 1)./(-2));
    HT  = (HD .+ HV);
    #
    cin = 0; cout = 0; cutPointsf = Vector{Float64}();
    for i in 1:length(HT)
        if HT[i] > 1
            if i <= n_points
                cin = cin + 1;
            else
                cout = cout + 1;
                push!(cutPointsf, i);
            end
        end
    end
    lenghtSlicesf = diff(cutPointsf) .- 1;
    #
    TerminatorToOSCAR(HD, HV, HT, cin, cout, cutPointsf, lenghtSlicesf)
end

struct Objects3D
    preok
    infoDF
    inis
    noInits
end
function ObjSplitter3D(file, medZ, minZ=0.5*medZ, maxZ=1.5*medZ, nCh=1)
    #
    data2d, data, sizesPlanes= data2d_chunking(file, 3, 4, 2, 1, 5, 8, 9, 19, nCh); data2d_copy = deepcopy(data2d);
    #
    indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5; indexA = 6; indexB = 7; indexType = 8; indexTag = 8; indexData = 9; indexInZ = 10; indexInRawData = 11;
    #
    preok = Vector{}();
    cutPointsF = Vector{}(); coh = Vector{}();
    val = Vector{}();
    #
    borderEffectLow = 1 #minZ - 1;
    inis = 0
    noInits = zeros(1,3)
    counter = 0
    counter_cin = 0
    count_pearson = 0
    for i in 1:(size(data2d_copy,1)) - borderEffectLow
        sliceZ = i;
        for j in 1:size(data2d_copy[i],1)
            indexR = j;
            goAns = 0 #initiators
            # mean_diameter = sqrt(mean([mean(data2d_copy[i][:,indexArea]) for i in 1:length(data2d_copy)])/pi)*2
            disappear_init = checkIntit(data2d_copy,i,j)
            if any(disappear_init .!= 0)
                noInits = vcat(noInits, disappear_init)
            end
            # println(noInits)
            if data2d_copy[i][j,indexTag] == 0
                inis += 1
                #
                preObj = preObjElongation(sliceZ,indexR, data2d_copy, maxZ, minZ, medZ);
                tempEl = preObj.tempEl; preVol = preObj.vVolV;
                # println("a'' ", tempEl[:,12])
                # println("b'' ", tempEl[:,13])
                distV = preObj.vDist; angV = preObj.vAng;
                pearson = preObj.vPearson;
                #
                cutPoint = size(tempEl, 1); factor = (pearson);
                if size(tempEl, 1) >= Int(round(minZ)) #&& factor > 0
                    #
                    if size(tempEl, 1) > Int(round(medZ))
                        # Esta es la version original

                        infoTerm = TerminatorReturns(distV[2:end], angV, medZ);
                        cin = infoTerm.numbIn; cout = infoTerm.numbOut;
                        cin = cin / medZ; cout = cout / (size(tempEl,1) - medZ + 1);
                        cutPV = infoTerm.cutOut; lengthPV = infoTerm.lengthOut;


                        if cin < cout && cout > 0 #cout > 0
                            counter_cin += 1
                            n_pointH = Int(cutPV[1] + 1); #H has one less element than var
                            #
                            if n_pointH -2 < 0
                                pearsonIn = 0; #corta por lo que definas antes del cutPoint
                            else
                                pearsonIn = getCovariance(distV[2:n_pointH-1], angV[1:n_pointH-1-1]);
                            end
                            #
                            if n_pointH -1 < 0 #Esto nunca deberia pasar.
                                pearsonOut = 0; #corta por lo que definas antes del cutPoint
                            else
                                pearsonOut = getCovariance(distV[n_pointH:end], angV[n_pointH-1:end]);
                            end
                            #
                            if pearsonIn > pearson
                                count_pearson += 1
                                cutPoint = Int(n_pointH);
                                counter += 1
                            end
                            #
                        end

                        # Aquí empiezan las funciones de prueba
                        # cutPoint, doubtedPoint = Termination_Chow(distV,angV)
                        # if cutPoint == 0
                        #     pearsonTot = getCovariance(distV[2:end], angV)
                        #     pearsonTilCut = getCovariance(distV[2:doubtedPoint+1], angV[1:doubtedPoint])
                        #     if pearsonTilCut > pearsonTot
                        #         cutPoint = doubtedPoint
                        #     end
                        # end
                    end
                    #
                    push!(cutPointsF, cutPoint);
                    push!(val, size(tempEl, 1));
                    push!(coh, factor);
                    for zz in 1:cutPoint
                        data2d_copy[Int(tempEl[zz,indexInZ])][Int(tempEl[zz,indexData]),indexTag] += 1
                        tempEl[zz,indexTag] = data2d_copy[Int(tempEl[zz,indexInZ])][Int(tempEl[zz,indexData]),indexTag]
                    end
                    cutTempEl = tempEl[1:cutPoint,:]
                    push!(preok, cutTempEl)
                    #
                end
            end
        end
    end
    k = DataFrame(V = val, C = cutPointsF, pearson = coh)
    #
    return Objects3D(preok, k, inis, noInits), counter_cin, count_pearson
end

function Obj_3D_lateClassification(pathXY, pathObject, nCh, dirOut = "new--3D_OSCAR_"; preok=0, summary=1)
    if ispath(dirOut) == 0
        mkdir(dirOut);
    end
    #index para data2d
    indexChanelXYfile = 19
    # indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5;
    # indexA = 6; indexB = 7; indexTag = 8; indexData2d = 9; indexInZ = 10;
    # indexInRawData = 11;
    #index for Line: this is the column indexing for each variable
    indexType = 12;
    # indexX = 1; indexY = 2; indexZ = 3; indexArea = 4; indexAng = 5; indexA = 6; indexB = 7;
    # indxVol = 8; indexVx = 9; indexVy = 10; indexVz = 11;
    #index para objects3d
    indX = 1; indY = 2; indZ = 3; indAr = 4; indAng = 5; indA = 6; indB = 7;
    data2d_index = 8; Zindex_data2d = 9; rawData_index = 10;
    indexCh_in = 11; # this must be always be at least 1 (even you don't have channels)
    #
    obj_dir = readdir(pathObject);
    # xyfiles = readdir(pathXY);
    for i in 1:length(obj_dir)
        # check if the list elment is a 3dobject folder
        tF = obj_dir[i]; #println("a ", file[end-3:end])
        obj_dir_name = "$pathObject\\$tF"
        println(obj_dir_name)
        if tF[1:10] == "3Dobjects_" && isdir(obj_dir_name)
            #abro la carpeta
            if ispath("$dirOut\\"*obj_dir[i]) == false && preok == 1
                mkdir("$dirOut\\"*obj_dir[i]);
            end
            objFiles = readdir(obj_dir_name);
            nameF = tF[11:end]*".txt";
            nameF = "$pathXY\\$nameF"
            if isfile(nameF)
                xy_file = readdlm(nameF, '\t', header = true)[1]
                summFile  =  zeros(1, indexType+nCh-1+1)
                for j in 1:length(objFiles)
                    a_t = objFiles[j];
                    tempPath = obj_dir_name*"\\$a_t";
                    #abro el carchivo
                    obj_t = readdlm(tempPath, '\t', header = true)[1];
                    #cambio de los channels
                    for k in 1:size(obj_t)[1]
                        # zz_t = Int.(obj_t[i, Zindex_data2d]);
                        # indIn_t = Int.(obj_t[i, data2d_index]);
                        for ch in 1:nCh
                            # obj_t[k, indexCh_in+j-1] = data2d_copy[zz_t][indIn_t, indexType+j-1]
                            obj_t[k, Int(indexCh_in+ch-1)] = xy_file[Int(obj_t[k, Int(rawData_index)]),Int(indexChanelXYfile+ch-1)]
                        end
                    end
                    if preok==1
                        # names!(obj,[:xc, :yc, :zc, :area, :angle, :majorc, :minorc, :type, :indData2D, :indZdata2D, :indexData])
                        objS = DataFrame(obj_t);
                        CSV.write("$dirOut\\"*obj_dir[i]*"\\$a_t", objS, delim = '\t')
                    end
                    if summary==1
                        ht = objectToRow_3D(obj_t, nCh);
                        # println("y antes push! ", size(summFile)); println(size(ht));
                        summFile = vcat(summFile, ht);
                        # push!(summFile, ht);
                        # println("y push! ", size(summFile));
                    end
                    #
                end
                if summary==1
                    # objects3D_df = summary_lateClass(preok, nCh)
                    # println("ole");
                    outPut = summFile[2:end,:];
                    println("save it");
                    # outPut = DataFrame(summFile);
                    # println("type ", typeof(outPut)); println(size(outPut))
                    objects3D_df = convert(DataFrame, outPut)
                    names = [:Xcenter, :Ycenter, :Zcenter, :AngleXY, :A, :B, :C, :Vol, :Vx, :Vy, :Vz]
                    for j in 1:nCh
                        names = vcat(names, Symbol("Channel$j"))
                    end
                    # names = vcat(names, Symbol("Type"))
                    names!(objects3D_df, names)
                    CSV.write("$dirOut\\Summary_"*tF[11:end]*".txt", objects3D_df, delim = '\t')
                end
            end
        end
    end
end

function sortPoints_ascending_mapToZ(xref,yref,xar,yar,zs)
    dist = zeros(length(xar))
    for i in 1:length(xar)
         dist[i] = sqrt((xar[i].-xref).^2 .+ (yar[i].-yref).^2)
        #dist[i] = ((xar[i].-xref).^2 .+ (yar[i].-yref).^2)
    end
    indicesInEachZ = vcat(collect(1:findfirst(zs .!= zs[1])-1),collect(1:length(zs)-(findfirst(zs .!= zs[1])-1)))
    inMinSorted = sortperm(vec(dist))
    indicesInEachZ = indicesInEachZ[inMinSorted]
    B = dist[inMinSorted]
    Zord = Int.(zs[inMinSorted])
    #dist = B
    return inMinSorted, B, dist, Zord, indicesInEachZ
end

function checkIntit(data2d_copy,i, j)
    indexX = 1; indexY = 2; indexArea = 3; indexAng = 4; indexZ = 5; indexA = 6; indexB = 7; indexType = 8; indexTag = 8; indexData = 9; indexInZ = 10; indexInRawData = 11
    # mean_diameter = sqrt(mean([mean(data2d_copy[i][:,indexArea]) for i in 1:length(data2d_copy)])/pi)*2
    # min_diameter = sqrt(mean([findmin(data2d_copy[i][:,indexArea])[1] for i in 1:length(data2d_copy)])/pi)*2
    # mean_diameter = sqrt(findmin(data2d_copy[i][:,indexArea])[1]/pi)*2
    # mean_diameter = sqrt(mean(data2d_copy[i][:,indexArea])/pi)
    mean_diameter = 0
    xc = data2d_copy[i][j,indexX]
    yc = data2d_copy[i][j,indexY]
    if i != 1
        xar = vcat(data2d_copy[i][:,indexX],data2d_copy[i-1][:,indexX])
        yar = vcat(data2d_copy[i][:,indexY],data2d_copy[i-1][:,indexY])
        zs = vcat(data2d_copy[i][:,indexInZ],data2d_copy[i-1][:,indexInZ])
        inMinSorted, B, dist, Zord, indicesInEachZ = sortPoints_ascending_mapToZ(xc,yc,xar,yar,zs)
        inMinSorted = [inMinSorted[k] for k in 1:length(inMinSorted) if data2d_copy[Zord[k]][indicesInEachZ[k],indexTag] == 1 && B[k] != 0]
        B = dist[inMinSorted]
        if length(B) != 0 && any(B .< mean_diameter)
            data2d_copy[i][j,indexTag] = 1
            disappear_init = hcat(data2d_copy[i][j,indexX], data2d_copy[i][j,indexY], data2d_copy[i][j,indexZ])
        end
    else
        xar = data2d_copy[i][:,indexX]
        yar = data2d_copy[i][:,indexY]
        inMinSorted, B, dist = sortPoints_ascending(xc,yc,xar,yar)
        # inMinSorted = inMinSorted[2:end]
        inMinSorted = [inMinSorted[k] for k in 1:length(inMinSorted) if data2d_copy[i][inMinSorted[k],indexTag] == 1 && B[k] != 0]
        B = dist[inMinSorted]
        if length(B) != 0 && any(B .< mean_diameter)
            data2d_copy[i][j,indexTag] += 1
            disappear_init = hcat(data2d_copy[i][j,indexX], data2d_copy[i][j,indexY], data2d_copy[i][j,indexZ])
        end
    end
    try
        return disappear_init
    catch
        disappear_init = zeros(1,3)
        return disappear_init
    end
end
