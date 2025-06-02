setBatchMode("hide"); 
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack limit redirect=None decimal=2");
//
setForegroundColor(255, 255, 255);
input = getDirectory("Input...");
files = getFileList(input);
ch = getNumber("Enter the nuclei channel number:", 1);
outputXY = input+"xy_files"+File.separator;  
File.makeDirectory(outputXY);
outputXZ = input+"xz_files"+File.separator; 
File.makeDirectory(outputXZ);
for (i = 0; i < files.length; i++) {
	a = endsWith(files[i], ".tif");
	if (a == 1){
    	run("Close All"); run("Clear Results");
    	print("Processing file index: " +files[i]);
    	file_path = ""+input+"/"+files[i];
		open(""+file_path+"");
		width = height = channels = slices = frames = -1; 
		getDimensions(width, height, channels, slices, frames);
		getDimensions(width, height, channels, slices, frames); 
		name1=File.nameWithoutExtension; rename(name1);
		if (bitDepth() != 8){
    		run("8-bit");
		}
		if (is("global scale") == 0){ //puse global scale = 0
    		selectWindow(""+name1+"");
    		run("Set Scale...", "distance=0 known=0 unit=pixel"); //quite pixel=1
    		selectWindow(""+name1+""); 
			getVoxelSize(dx, dy, dz, unit);
		}else{	
			selectWindow(""+name1+""); 
			getVoxelSize(dx, dy, dz, unit);
			}
		if (channels > 1) {
			selectWindow(""+name1+""); run("Split Channels"); 
			selectWindow("C"+ch+"-"+name1+""); rename("nuclei");
		}else{
			selectWindow(""+name1+"");  rename("nuclei");
		}
		selectWindow("nuclei"); setThreshold(1,255); selectWindow("nuclei");
		setOption("BlackBackground", true); run("Convert to Mask", "method=Default background=Dark black");
		selectWindow("nuclei"); run("Watershed", "stack");
		selectWindow("nuclei"); run("Analyze Particles...", "size=0-Infinity display stack");
		selectWindow("nuclei"); resetThreshold();
	
		slice = newArray(nResults); xc = newArray(nResults); yc = newArray(nResults);
		area = newArray(nResults); angle = newArray(nResults); type = newArray(nResults); 
		major = newArray(nResults);  minor = newArray(nResults);
		intDen = newArray(nResults); PromInt = newArray(nResults); 
		feret = newArray(nResults); xfer = newArray(nResults); yfer = newArray(nResults); 
		minferet = newArray(nResults); anglefer = newArray(nResults); 
		xb = newArray(nResults); yb = newArray(nResults); wb = newArray(nResults); hb = newArray(nResults);
		areaNuclei = newArray(nResults);
		inA=0; newImage("base2", "8-bit grayscale-mode", width, height, 1, slices, 1);
		for(r=0;r<nResults;r++){
			aell = getResult("Area",r);	
			rad = getResult("Minor",r); xref = getResult("X",r);
			yref = getResult("Y",r); zref = getResult("Slice",r);
			a = getResult("Major", r)/2; b = getResult("Minor", r)/2;
			angleXY = getResult("Angle", r); angref=angleXY*(PI/180)*(-1);
			xv1 = xref + (cos(angref) * a); yv1 = yref + (sin(angref) * a);
			xv2 = xref - (cos(angref) * a); yv2 = yref - (sin(angref) * a);					
			
			selectWindow("nuclei"); makeEllipse(xv1, yv1, xv2, yv2, b/a);
			selectWindow("nuclei"); setThreshold(1, 255);
			aell=getValue("Area limit");
			selectWindow("nuclei"); resetThreshold();
			selectWindow("nuclei"); run("Select None");
			no = aell;
			
			if (channels > 1) {
				channel_array = newArray(channels);
				for (i = 0; i < channels; i++) {
    				channel_array[i] = i + 1;
				}
				newChannels = newArray();
				for (i = 0; i < lengthOf(channel_array); i++) {
    				if (channel_array[i] != ch) {
        				newChannels = Array.concat(newChannels,channel_array[i]);
    				}
				}
				n_channels = newArray();
				//keys = newArray(newChannels);
				for (i = 0; i < newChannels.length; i++){
					key = "n_" + newChannels[i];
					selectWindow("C"+newChannels[i]+"-"+name1+""); setSlice(zref);
					selectWindow("C"+newChannels[i]+"-"+name1+""); makeEllipse(xv1, yv1, xv2, yv2, b/a);
					selectWindow("C"+newChannels[i]+"-"+name1+""); n_channels[key]=getValue("RawIntDen");
				}
			}
			selectWindow("base2"); setSlice(zref);
			selectWindow("base2"); makeEllipse(xv1, yv1, xv2, yv2, b/a);
			selectWindow("base2"); run("Fill", "slice");
			selectWindow("base2"); run("Select None");	
			
			slice[inA]=getResult("Slice",r); 
			area[inA]=getResult("Area",r);
			xc[inA]=getResult("X",r);
			yc[inA]=getResult("Y",r);
			major[inA]=getResult("Major",r); 
			minor[inA]=getResult("Minor",r); 
			angle[inA]=getResult("Angle",r);
			intDen[inA]=getResult("IntDen",r);
			PromInt[inA]=getResult("Mean",r);
			feret[inA]=getResult("Feret",r); 
			xfer[inA]=getResult("FeretX",r); 
			yfer[inA]=getResult("FeretY",r); 
			anglefer[inA]=getResult("FeretAngle",r); 
			minferet[inA]=getResult("MinFeret",r);
			xb[inA]=getResult("BX",r); 
			yb[inA]=getResult("BY",r); 
			wb[inA]=getResult("Width",r); 
			hb[inA]=getResult("Height",r);
			areaNuclei[inA]=no;
			if (channels >1) {
				intensityValues = newArray(lengthOf(newChannels));
				for (i = 0; i < lengthOf(newChannels); i++) {
    				key = "n_" + newChannels[i];
    				intensityValues[i] = n_channels[key];  // Retrieve stored intensity
				}

				// Add the intensity values dynamically to the table
			for (i = 0; i < lengthOf(newChannels); i++) {
   				channelKey = "intCh" + newChannels[i]; // e.g., intCh1, intCh2, etc.
    			setResult(channelKey, inA, intensityValues[i]);
				}
			}
			
			inA=inA+1;
						
		}
		run("Clear Results");
		Table.create("2Dinf");
		Table.setColumn("Slice", slice);
		Table.setColumn("Area", area);
		Table.setColumn("xc", xc);
		Table.setColumn("yc", yc);				
		Table.setColumn("angle (degrees)",angle);
		Table.setColumn("intDen", intDen);
		Table.setColumn("PromInt", PromInt);
		Table.setColumn("majorc", major);
		Table.setColumn("minorc", minor);		
		Table.setColumn("Feret", feret);
		Table.setColumn("FeretX", xfer);
		Table.setColumn("FeretY", yfer);
		Table.setColumn("FeretAngle", anglefer); 
		Table.setColumn("MinFeret", minferet);
		Table.setColumn("BX", xb);
		Table.setColumn("BY", yb);
		Table.setColumn("BWidth", wb); 
		Table.setColumn("BHeight", hb);
		Table.setColumn("Area_Nulcei", areaNuclei);
		// Dynamically generate Intensity_ChannelN columns
		if (channels > 1) {

			for (i = 0; i < lengthOf(newChannels); i++) {
	    		channelNumber = newChannels[i]; // Get channel number (e.g., 1, 2, 3)
	    		columnName = "Intensity_Channel" + channelNumber; // e.g., "Intensity_Channel1"
	  	    	// Create an array to store intensity values for this channel
	   			intensityValues = newArray(inA);
	 	   		for (j = 0; j < inA; j++) {
	        		key = "n_" + channelNumber;
	        		intensityValues[j] = n_channels[key]; // Retrieve stored intensity
	    		}
	    		// Add the intensity values to the table
	    		Table.setColumn(columnName, intensityValues);
			}
		}

		nameFileR=""+name1+"";
		Table.save(outputXY + nameFileR + ".txt");
		Table.reset("2Dinf");
			
		minRz = 3;
		selectWindow("base2");getVoxelSize(dx, dy, dz, unit);
		selectWindow("base2"); run("Reslice [/]...", "output="+dz+" start=Top avoid");
		selectWindow("Reslice of base2"); rename("xz_base");
		selectWindow("xz_base"); setThreshold(1,255); 
		setOption("BlackBackground", true); selectWindow("xz_base");
		run("Convert to Mask", "method=Default background=Dark black");
		selectWindow("xz_base"); run("Minimum...", "radius="+minRz+" stack");
		selectWindow("xz_base"); run("Maximum...", "radius="+minRz+" stack");						
		selectWindow("xz_base"); setThreshold(1,255); 
		setOption("BlackBackground", true); selectWindow("xz_base");
		run("Convert to Mask", "method=Default background=Dark black");
		selectWindow("xz_base"); run("Watershed", "stack");
		selectWindow("xz_base"); setThreshold(1,255); 
		run("Clear Results");
		run("Analyze Particles...", "size=0-Infinity display stack");
		selectWindow("xz_base"); resetThreshold();
		sliceZ = newArray(nResults); xcZ = newArray(nResults); ycZ = newArray(nResults); 
		areaZ = newArray(nResults); angleZ = newArray(nResults); 
		majorZ = newArray(nResults);  minorZ = newArray(nResults); 
		feretZ = newArray(nResults); xferZ = newArray(nResults); yferZ = newArray(nResults); 
		minferetZ = newArray(nResults); angleferZ = newArray(nResults); 
		xbZ = newArray(nResults); ybZ = newArray(nResults); wbZ = newArray(nResults); hbZ = newArray(nResults);
		inA=0;
		for(r=0;r<nResults;r++){
			sliceZ[inA]=getResult("Slice",r); 
			areaZ[inA]=getResult("Area",r);
			xcZ[inA]=getResult("X",r);
			ycZ[inA]=getResult("Y",r);
			majorZ[inA]=getResult("Major",r); 
			minorZ[inA]=getResult("Minor",r); 
			angleZ[inA]=getResult("Angle",r);
			feretZ[inA]=getResult("Feret",r); 
			xferZ[inA]=getResult("FeretX",r); 
			yferZ[inA]=getResult("FeretY",r); 
			angleferZ[inA]=getResult("FeretAngle",r); 
			minferetZ[inA]=getResult("MinFeret",r);
			xbZ[inA]=getResult("BX",r); 
			ybZ[inA]=getResult("BY",r); 
			wbZ[inA]=getResult("Width",r); 
			hbZ[inA]=getResult("Height",r);
			inA=inA+1;
		}
		run("Clear Results");
		Table.create("2Dinf_XZ");
		Table.setColumn("Slice", sliceZ);
		Table.setColumn("Area", areaZ);
		Table.setColumn("xc", xcZ);
		Table.setColumn("yc", ycZ);				
		Table.setColumn("angle (degrees)",angleZ);
		Table.setColumn("majorc", majorZ);
		Table.setColumn("minorc", minorZ);
		Table.setColumn("Feret", feretZ);
		Table.setColumn("FeretX", xferZ);
		Table.setColumn("FeretY", yferZ);
		Table.setColumn("FeretAngle", angleferZ); 
		Table.setColumn("MinFeret", minferetZ);
		Table.setColumn("BX", xbZ);
		Table.setColumn("BY", ybZ);
		Table.setColumn("BWidth", wbZ); 
		Table.setColumn("BHeight", hbZ);
		nameFileR=""+name1+"";
		Table.save(outputXZ + nameFileR + ".txt");
		Table.reset("2Dinf_XZ");
		selectWindow("base2"); close();
		nameIm=""+name1+""; 	
		selectWindow("nuclei");	
		close(); 		
	}
}
//setBatchMode("exit and display");
