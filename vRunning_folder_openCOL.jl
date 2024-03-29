# include(pathMain*"\\SysBio3D_methods.jl")
#widgets: que hubiera algo rollo input(); en fiji, en jupyter notebook, si no cambiar a mano...
# path = "C:\\Users\\i7\\Desktop\\TestFinal\\input\\txtFiles";
# path = "C:\\Users\\i7\\Desktop\\Nueva carpeta33";
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\Outputs\\noSegmented"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\Outputs\\Watersheed"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\Outputs\\vHomeMade"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\Outputs\\watersheed3d"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\Outputs\\txt"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\n\\txt"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc_reales\\Input_images\\Processed\\2dsegm"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc_reales\\Input_images\\processed_varRads\\2dsegm"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc_reales\\Input_images\\processed_varRads\\3dsegm"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc_1500_n-3"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\n-1_n-2\\n-2"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\n\\txt"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\noiseArounXYZ\\v2\\corwded_processing";
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\noiseArounXYZ\\v2\\maxmin_processing";
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\noiseArounXYZ\\v2\\max3\\blurreo_processing\\segm3d";
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\noiseArounXYZ\\v2\\max1\\noCrowded\\segm2d"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\nc\\noiseArounXYZ\\v2\\max1\\noCrowded\\segm3d";
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\PAPER_OSCAR\\fig2\\new-Output\\rad2.5"
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\PAPER_OSCAR\\fig2\\recort\\raw\\fr\\nowaterInz"
include("C:\\Users\\i7\\Desktop\\working_OSCAR\\vPearson_outliers_openCOL.jl")
# path = "C:\\Users\\mario\\OneDrive\\Escritorio\\PAPER_OSCAR\\fig2\\ultima"
# path = "C:\\Users\\i7\\Desktop\\PruebaProcSpin\\imagen_processings\\processed_images_36_2"
# path = "G:\\Mi unidad\\Proyectos\\colaboracionKrotenberMario\\processing\\dapi"
#path = "G:\\Mi unidad\\Diego\\Drogas\\SelectedProcessings\\Cyc\\processed"
path = "G:\\Mi unidad\\Diego\\Drogas\\ultimas\\xav_selected\\processed"
# dirspath = "G:\\Mi unidad\\Proyectos\\Oscar\\drogas_2023_todas"
# \\outcomes
dirspath = "G:\\Mi unidad\\Proyectos\\Oscar\\drogas_2023_todas\\outcomes\\DAPT\\processing"
dirspath = "G:\\Mi unidad\\Proyectos\\Oscar\\drogas_2023_todas\\SU_greyscale_compr\\Processings_AutoMinR\\3D_OSCAR_objects\\processing\\Processings_AutoMinR\\3D_OSCAR_objects\\processing"
dirspath = "G:\\Mi unidad\\Proyectos\\Oscar\\AQUI ESTA LO DE LA TESIS PATENTE PAPER (2020-21)\\images-hom\\images"

# dirspath = "G:\\Mi unidad\\Diego\\29_12_23\\recortadas"
dirs = readdir(dirspath)
for i in 86:86
    path = dirspath*"\\"*dirs[i]
    # "\\Processings_AutoMinR\\"
    nameFolO = "\\3D_OSCAR_objects"
    a = Date(now()); fileOutput = path*nameFolO;
    if ispath(fileOutput) == 0
        mkdir(fileOutput);
        #rm(fileOutput, recursive = true)
    end
    inputPath = path*"\\xy_files";
    files = readdir(inputPath)
    #
    names2 = Vector{String}(); numb = Vector{}(); zminv = Vector{}();
    zmedv = Vector{}(); zmaxv = Vector{}(); inisv = Vector{}(); times = Vector{}()
    rad = [6,9,12]; nCh = 3; #
    # for i in length(files):-1:1
    for i in 1:1
        # stats = @timed begin
        # @time begin
        file  = files[i]; println("name ", file)
        fileXZ = filter(fileXZ->occursin(fileXZ,file), readdir(path*"\\xz_files"));
        namez = string(fileXZ[1])
        # namez = string(fileXZ[1])*"_xz"
        pathXZ = path*"\\xz_files\\"*string(fileXZ[1]);
        zp = min_max_Zdef(pathXZ);
        # zp = Zmin_Zmed_Zmax2(pathXZ);
        # minZinterval = zp.minZ; maxZinterval = zp.maxZ; medZinterval = zp.medZ;
        minZinterval = floor(zp.minZ/2); maxZinterval = floor(zp.minZ+zp.minZ/2); medZinterval = zp.minZ; #Esta línea es la normal
        # maxZinterval = zp.minZ; medZinterval = round(Int,zp.minZ/2); minZinterval = round(Int,zp.minZ/4) #Para imagenes con poco slices
        # medZinterval = 30; minZinterval = floor(medZinterval/2); maxZinterval = floor(medZinterval+medZinterval/2); #Solo para la mosca
        println("zp is: ", (minZinterval,medZinterval,maxZinterval))
        # minZinterval = 3; maxZinterval = 9; medZinterval = 5;
        # println("zp is: ", zp)
        # minZinterval = 1; medZinterval = 5; maxZinterval = 9;
        if minZinterval < 1
            println("z param tends to 1. This value means that all 2D objects are 3D objects.")
            push!(names2, file); push!(numb, 000); push!(zminv, 000); push!(zmedv, 000); push!(zmaxv, 000);
        else
            @time begin
                folder = inputPath*"\\$file";
               global objc3d, counter_cin, counter_pearson = ObjSplitter3D(folder, medZinterval, minZinterval, maxZinterval, nCh);
            end
            # ft = file[1:end-4]
            # output_obj = fileOutput*"\\3Dobjects_$ft";
            # if ispath(output_obj) == 0
            #     mkdir(output_obj);
            # end
            # preoObj_df= toDataFrame(objc3d, output_obj, file, nCh);

            h = objc3d.infoDF
            CSV.write("$fileOutput\\info_$file.txt", h, delim = '\t')
            # #
            objects = objc3d.preok;
            println("hay ", size(objects,1))
            summary_df = summary_objs(objc3d, nCh); #summary_df = hcat(summary_df, h);
            # println(summary_df)
            CSV.write("$fileOutput\\Summary_$file", summary_df, delim = '\t')

            push!(names2, file);
            push!(numb, size(objc3d.preok, 1));
            push!(zminv, minZinterval);
            push!(zmedv, medZinterval);
            push!(zmaxv, maxZinterval);
            push!(inisv, objc3d.inis);
        end
    end
    # push!(times, stats.time)
    end
    col_labels = ["name", "obj", "zmin", "zmed", "zmax", "inis", "computation_time"];
    col_labels = Symbol.(col_labels);
    df = DataFrame(hcat(names2, numb, zminv, zmedv, zmaxv, inisv, times), col_labels);
    CSV.write("$fileOutput\\BasicInfo_Image(s).txt", df, delim = '\t')
end
