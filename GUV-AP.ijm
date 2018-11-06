/*
+++++++++++++++++++++++++++++++++++++++++++++++++
+  	Automated GUV analysis software 	+
+	Develloped by Taras Sych					+
+ 	This plug-in is designed to investigate 	+
+	membrane dye partitioning and protein 		+
+ 	binding either to membrane lipids and/or	+
+	proteins. In addition, protein binding is 	+
+	evaluate based on the fluorescence			+
+	intensity method. 							+
+   Version  --- 2017							+
+												+	
+++++++++++++++++++++++++++++++++++++++++++++++++

+
Image recommandation
•	Raw data gray scale image (no scale or color bars)
•	The number of channels 1 - 3
•	Image can contain several frames/slices, but only first frame and/or slice will be processed.

*/

requires("1.51p");
ScreenClean();
getDateAndTime(year, month, dayofWeek, dayofMonth, hour, minute, second, msec);


//---------------------------------------Input parameter section--------------------------------------------
//----------------------------------------------------------------------------------------------------------
	tc = 0.05;				// threshold coeffitient for GUV detection

	min_detected =200; 		// minimal area of detectede particles in pixels

	ce = 1.2; 				// circle extension

	cirularity_check = 0.95;	// minimal possible circularity
	
	file_extension = "tif"; // extension of processed files (NO DOT!)


	alpha_step = 6; 		// sampling angle step

	setBatchMode(false);	// activate/deactivate BatchMode
//----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------


//----------------Select directory containing data to process-----------------------------------------------
//----------------------------------------------------------------------------------------------------------
dir = getDirectory("Choose directory");
list = getFileList(dir);
//----------------------------------------------------------------------------------------------------------
ttt=tc;


//----------------Create directory for data output----------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
dir1 = dir + " " + year +"." + month + "." + dayofMonth + " " + hour + "." + minute + " - Result\\" ;
File.makeDirectory(dir1);
//----------------------------------------------------------------------------------------------------------

title0 = "membrane";
title1 = "protein_1";
title2 = "protein_2";


//----------------Select channel to process (GUI)-----------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
Dialog.create("Analyze channel");
Dialog.addCheckbox("channel 1", true);
Dialog.addString("Name :", title0);
Dialog.addCheckbox("channel 2", false);
Dialog.addString("Name :", title1);
Dialog.addCheckbox("channel 3", true);
Dialog.addString("Name :", title2);
Dialog.show();

ch1 = Dialog.getCheckbox();
title0 = Dialog.getString();

ch2 = Dialog.getCheckbox();
title1 = Dialog.getString();

ch3 = Dialog.getCheckbox();
title2 = Dialog.getString();

i = 0;
titles = newArray (3);

if (ch1== true) {
	titles [i] = title0;
	i++;
}
if (ch2 == true) {
	titles [i] = title1;
	i++;
}
if (ch3 == true) {
	titles [i] = title2;
	i++;
}
channel_counter=i;

background = newArray (channel_counter);
background_noise_arr = newArray (channel_counter);
//----------------------------------------------------------------------------------------------------------

//----------------Select channel for GUV detection and correction method------------------------------------
//----------------------------------------------------------------------------------------------------------
		
		  html = "<html><font size=-1><b>Circular check:</b> Non-phase separated GUVs<br> <b>GUV stitching:</b>  Phase separated GUVs <br></font>";
		Dialog.create("GUV detection");
		 Dialog.addChoice("Channel used for Segmentation:", newArray(title0, title1, title2,"Sum all channels"));
	     	 Dialog.addRadioButtonGroup("Segmentation mode",newArray("Circular check (Non-phase separated)","GUV stitching (Phase separated)"), 2, 1,"Circular check (Non-phase separated)");
	     	 Dialog.addHelp(html);
	    Dialog.addChoice("Preview:", newArray("None","First file", "Select file"));
		Dialog.show();
	
		channel = Dialog.getChoice();
		 if (channel=="Sum all channels")
		   channel="Overlay";
	 	SegMode= Dialog.getRadioButton();
	 	preview_mode = Dialog.getChoice();

if (preview_mode == "First file" || preview_mode == "Select file") {
setBatchMode(false);
}


if (SegMode == "Circular check (Non-phase separated)") {
	circ_check = true;
	GUV_stitch = false;
}

if (SegMode == "GUV stitching (Phase separated)") {
	circ_check = false;
	GUV_stitch = true;
}
//----------------------------------------------------------------------------------------------------------

//----------------Create Parameter_log.txt and save parameters----------------------------------------------
//----------------------------------------------------------------------------------------------------------
xl2=File.open(dir1+ File.separator + "Parameter_log.txt");

print (xl2, "ImageJ version: " + getVersion());
print (xl2," ");
print (xl2, "\t\t" + "GUV detection parameters");
print (xl2," ");
print (xl2, "Channel for ROI definition: " + channel);
print (xl2, "Thresholding for ROIs:" + " tc = " + ttt);
print (xl2, "Method for GUV detection: " + "Fit_Circle");
if (circ_check == true) print (xl2, "Circular check");
if (GUV_stitch == true) print (xl2, "GUV stitching");
print (xl2, "Circle extension coefficient: " + "ce = " + ce);
print (xl2," ");
print (xl2, "\t\t" + "Calculation paramters");
print (xl2," ");
print (xl2, "Angle step: alpha_step = " + alpha_step + " deg");
print (xl2," ");


File.close(xl2);
//----------------------------------------------------------------------------------------------------------

//----------------Various arrays initialization-------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
Total_membrane = newArray (361);
Total_protein_1 = newArray (361);
Total_protein_2 = newArray (361);

Max_intensity_membrane = newArray (360/alpha_step + 1);
Max_intensity_protein_1 = newArray (360/alpha_step + 1);
Max_intensity_protein_2 = newArray (360/alpha_step + 1);

Angle_array = newArray (361);

counter = 0;

summary_ch1 = newArray;

summary_ch2 = newArray;

summary_ch3 = newArray;
//----------------------------------------------------------------------------------------------------------

//----------------Check data folder for suitabke files------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
file_counter = 0;
for (number_of_file = 0; number_of_file<list.length; number_of_file++){
	if (endsWith(list[number_of_file],file_extension) == 1){
		file_counter ++;
	}
}
//----------------------------------------------------------------------------------------------------------

if (file_counter == 0){
	waitForUser("Error", "No files detected. Plese, check extension");
} else {
//----------------Loop analyses files one by one------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------	
	if (preview_mode == "Select file") {
		file = File.openDialog("Choose a File");
		run("Bio-Formats Importer", "open=file autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		image123 = getTitle();
	}
	for (number_of_file = 0; number_of_file<list.length; number_of_file++){
		if (endsWith(list[number_of_file],file_extension) == 1){

		//----------File import


		
		file = dir + list[number_of_file];
		//file_marked = dir1 + File.separator + list[number_of_file];

		if (preview_mode == "None" || preview_mode == "First file"){
		run("Bio-Formats Importer", "open=file autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		} else selectWindow (image123);
		
		file_marked = dir1 + File.separator + File.nameWithoutExtension;
		file_excell = dir1 + File.separator + File.nameWithoutExtension;
		run("32-bit");

		rename("raw_data_0");
		run("Duplicate...", "title=raw_data duplicate frames=1");

		//----------File split on color channels
		if (channel_counter > 1) run("Split Channels");
		
		if (channel_counter > 1) selectWindow("C1-raw_data");
		rename(titles[0]);
		name_back = titles[0] + "_background";
		run("Duplicate...", " ");
		rename(name_back);
		
		if (channel_counter == 2){
			selectWindow("C2-raw_data");
			rename(titles[1]);
			name_back = titles[1] + "_background";
			run("Duplicate...", " ");
			rename(name_back);
			imageCalculator("Add create", titles[0],titles[1]);
			
			rename ("Overlay");
		}

		

		if (channel_counter == 3){
			selectWindow("C2-raw_data");
			rename(titles[1]);
			name_back = titles[1] + "_background";
			run("Duplicate...", " ");
			rename(name_back);
			selectWindow("C3-raw_data");
			rename(titles[2]);
			name_back = titles[2] + "_background";
			run("Duplicate...", " ");
			rename(name_back);
			imageCalculator("Add create", titles[0],titles[1]);
			
			rename ("Overlay");
			imageCalculator("Add create", "Overlay",titles[2]);
			selectWindow("Result of Overlay");
			rename ("Overlay");
		}

		//------Channel for ROI definition--------------------
		name = channel;

		Stack.getDimensions(width33, height33, channels33, slices33, frames33);

		getPixelSize(unit, pixelWidth, pixelHeight); 
		
		selectWindow(name);

		
		
		

		//-----Threshold applying and particle detection
		run("Duplicate...", " ");
		rename("mask");

		run("Despeckle");
		
		run("Set Measurements...", "min redirect=None decimal=3");
		run("Measure");

		top = getResult("Max",0);
		bottom = (top - getResult("Min",0)) * tc + getResult("Min",0);
		setAutoThreshold("Default dark");
		setThreshold(bottom, top);

		run("Set Measurements...", "area centroid bounding shape redirect=None decimal=3");
		run("Analyze Particles...", "size=min_detected-Infinity pixel show=Masks display exclude clear  add in_situ");


		//------Creating background images for each color channel
		for (chan = 0; chan < channel_counter; chan ++){
			name_back = titles[chan] + "_background";
			selectWindow (name_back);
			
			
			for (roi_i = 0; roi_i < roiManager("count"); roi_i ++){
				roiManager("Select", roi_i);
				run("Clear", "slice");
			}
		}

		

		//Smoothing and skeletonisation
		selectWindow ("mask");
		run("Select None");
		for (i=0; i<10; i++)
			run("Median...","radius=2");
			/*run("Smooth");
		setOption("BlackBackground", false);
		run("Make Binary");*/
		run("Skeletonize");


		if (roiManager("count") > 0){
				roiManager("Delete");
			}
		run("Clear Results");

		run("Set Measurements...", "area centroid bounding shape redirect=None decimal=3");
		run("Analyze Particles...", "size=5-Infinity pixel show=Masks display exclude clear  add in_situ");
		GUVs_total = nResults;
		
		//-----------------Circularity check and filtering
		if (circ_check == true){
			for (number_of_GUV = 0; number_of_GUV < GUVs_total; number_of_GUV ++){
				round123 = getResult("Round",number_of_GUV);
				if ( round123 < cirularity_check){
					IJ.deleteRows(number_of_GUV,number_of_GUV) ;
					roiManager("Select", number_of_GUV);
					roiManager("Delete");
					GUVs_total = GUVs_total - 1;
					number_of_GUV = -1;
				}
			}
		}

		
		R_array = newArray ( GUVs_total );
		x_c_array = newArray ( GUVs_total );
		y_c_array = newArray ( GUVs_total );

		for (i=0; i<GUVs_total; i++){
			run("Select None");
			roiManager("Select", i);
			getSelectionCoordinates(x_points, y_points); 
			makeSelection("point",x_points, y_points);
			run("Fit Circle");
			run("Clear Results");


			run("Set Measurements...", "centroid bounding redirect=None decimal=3");
			run("Measure");
			x_c_array [i] = getResult("X", 0);
			y_c_array [i] = getResult("Y", 0);

			width = getResult("Width", 0);
			height = getResult("Height", 0);

			R_array[i] = (width+height)/4;
		}

		//------Phase separated GUV stitching
		kk = 0;
		stitched_array = newArray(1000);
		if (GUV_stitch == true){
			GUV_counter = 0;

		for (i=0; i<GUVs_total; i++){
			q=newArray(1);
			q[0]=i;
			stitched = false;
			roi_counter = 1;
			roiManager("Sort");
			for(j=0; j<GUVs_total; j++){
				
				if(sqrt((x_c_array [i] - x_c_array [j])*(x_c_array [i] - x_c_array [j]) + (y_c_array [i] - y_c_array [j])*(y_c_array [i] - y_c_array [j]))<(R_array[i]+R_array[j]) && j!=i){
					

					q = Array.concat (q,j);
					roi_counter++;
					stitched = true;
	
				}
			}

				

				stitched_array [kk] = stitched;

				kk++;
				
				roiManager("Sort");
				roiManager("Select", q);
				if (roi_counter > 1) roiManager("Combine");
				roiManager("Add");
				Num = roiManager("count");
				roiManager("Select", Num -1);
				roi_name = "GUV_" + GUV_counter;
				roiManager("Rename", roi_name);

			
				roiManager("Select", q);
				roiManager("Delete");

				roiManager("Sort");


				x_c_array = delete_array_elements (x_c_array,q);
				y_c_array = delete_array_elements (y_c_array,q);
				R_array = delete_array_elements (R_array,q);

				


				
				GUV_counter ++;

				
				
				GUVs_total = GUVs_total-roi_counter;
				i=-1;

		}


		R_array1 = newArray ( GUV_counter );
		x_c_array1 = newArray ( GUV_counter );
		y_c_array1 = newArray ( GUV_counter );
		
		for (i=0; i< GUV_counter; i++){
			
			run("Select None");
			roiManager("Select", i);


			if (stitched_array [i] == 1){
			points = get_all_points_in_a_ROI();

			x_points = get_x(points);

			y_points = get_y(points); }
			
			if (stitched_array [i] == 0) {getSelectionCoordinates(x_points, y_points);}
			
			makeSelection("point",x_points, y_points);
			run("Fit Circle");
			run("Clear Results");

			setForegroundColor(255, 255, 255);
			run("Draw", "slice");

			run("Set Measurements...", "centroid bounding redirect=None decimal=3");
			run("Measure");
			x_c_array1 [i] = getResult("X", 0);
			y_c_array1 [i] = getResult("Y", 0);

			width = getResult("Width", 0);
			height = getResult("Height", 0);

			R_array1[i] = (width+height)/4;
		}

				x_c_array = x_c_array1;
				y_c_array = y_c_array1;
				R_array = R_array1;

				

				GUVs_total = GUV_counter;
			
		} 

	

		//-----------Marked GUVs image saving

		file_marked = dir1 + File.separator + list[number_of_file];
		

		selectWindow ("raw_data_0");
		run("Duplicate...", "title=raw_data_saving duplicate frames=1");
		roiManager("Show None");
		roiManager("Show All");

		saveAs("Tiff", file_marked);
		close();

		
		//------Background calculation
		for (chan = 0; chan < channel_counter; chan ++){
			name_back = titles[chan];
			selectWindow(name_back);
			background_noise = 0;
			for(i=0; i<GUVs_total; i++){
				makeOval ((x_c_array[i] - R_array[i]/3)/pixelWidth, (y_c_array[i] - R_array[i]/3)/pixelWidth, 2*R_array[i]/3/pixelWidth, 2*R_array[i]/3/pixelWidth);
				
				run("Set Measurements...", "mean min redirect=None decimal=3");
				run("Measure");
				background_noise = background_noise + getResult("Mean",0);
				run("Clear Results");
				run("Select None");
			}
			
			background_noise = background_noise/GUVs_total;
			background_noise_arr [chan] = background_noise;

			
		}
		
		//------Cleaning background images and free protein calculation
		for (chan = 0; chan < channel_counter; chan ++){
			name_back = titles[chan] + "_background";
			selectWindow(name_back);

			for(i=0; i<GUVs_total; i++){
				makeOval ((x_c_array[i] - R_array[i]*ce)/pixelWidth, (y_c_array[i] - R_array[i]*ce)/pixelWidth, 2*R_array[i]*ce/pixelWidth, 2*R_array[i]*ce/pixelWidth);
				run("Clear", "slice");
				run("Select None");
			}
			
			run("Set Measurements...", "mean min redirect=None decimal=3");
			run("Measure");
			top = getResult("Max",0);
			run("Clear Results");
			//setAutoThreshold("Default dark");
			setThreshold(0.01, top);
			run("NaN Background", "stack");

			run("Set Measurements...", "mean min redirect=None decimal=3");
			run("Measure");
			
	
			background [chan] = getResult("Mean",0);
			if (background [chan] != background [chan] || background [chan] == 0) background [chan] = 1;
			run("Clear Results");

			selectWindow(name_back);
			
			close();
		
		}

		
	
		
		//--------------Preparation for processing-------------------
		//--------------channel 1-------------------
		chan = 0;
		if (ch1 == true){
		
			//-------Threshold
			selectWindow(title0);
			run("Set Measurements...", "mean min redirect=None decimal=3");
			run("Measure");
			
			top = getResult("Max",0);
			bottom = (background[chan] - background_noise_arr [chan])/2;

			bottom_ch1 = bottom;
			run("Clear Results");
			selectWindow(title0);
			setAutoThreshold("Default dark");
			//run("Threshold...");
			setThreshold(bottom, top);

			run("NaN Background", "stack");

			
		chan++;
		}
		
		//--------------channel 2-------------------
		
		if (ch2 == true){
			//-------Threshold
			selectWindow(title1);
			run("Set Measurements...", "mean min redirect=None decimal=3");
			run("Measure");

			top = getResult("Max",0);
			bottom = (background[chan] - background_noise_arr [chan])/2;

			bottom_ch2 = bottom;
			run("Clear Results");
			selectWindow(title1);
			setAutoThreshold("Default dark");
			//run("Threshold...");
			setThreshold(bottom, top);

			run("NaN Background", "stack");

			
		chan++;
		}
		
		//--------------channel 3-------------------

		if (ch3 == true){
			//-------Thresholding
			selectWindow(title2);
			run("Set Measurements...", "mean min redirect=None decimal=3");
			run("Measure");

			top = getResult("Max",0);
			bottom = (background[chan] - background_noise_arr [chan])/2;

			bottom_ch3 = bottom;
			run("Clear Results");
			selectWindow(title2);
			setAutoThreshold("Default dark");
			//run("Threshold...");
			setThreshold(bottom, top);

			run("NaN Background", "stack");

			
		chan++;
		}

//--------------loop for GUV processing------------------------------------------

		for (number_of_GUV = 0; number_of_GUV < GUVs_total; number_of_GUV ++){


			//Ring selection
				newImage("ring_mask", "8-bit black", width33, height33, 1);
				makeOval ((x_c_array[number_of_GUV] - R_array[number_of_GUV])/pixelWidth, (y_c_array[number_of_GUV] - R_array[number_of_GUV])/pixelWidth, 2*R_array[number_of_GUV]/pixelWidth, 2*R_array[number_of_GUV]/pixelWidth);
				
				setForegroundColor(255, 255, 255);
				
				run("Draw", "slice");
				run("Select None");
				run("Options...", "iterations=1 count=1");
				run("Erode");
				run("Erode");
				run("Erode");

				roiManager("Deselect");
				roiManager("Delete");
				run("Invert LUT");

				run("Analyze Particles...", "size=5-Infinity show=Masks exclude add");
				run("Invert LUT");
				run("Analyze Particles...", "size=5-Infinity show=Masks exclude add");
				close();
				close();
				roiManager("Select", newArray(0,1));
				close();
				roiManager("XOR");
				roiManager("Add");
				roiManager("Delete");
				roiManager("Select", 0);
				roiManager("Select", 0);
				roiManager("Rename", "Ring");
			//End of ring selection

			selectWindow("mask");
			run("Select None");
			/*roiManager("Deselect");
			roiManager("Delete");*/

			x_c = x_c_array [number_of_GUV];
			y_c = y_c_array [number_of_GUV];
			R = R_array [number_of_GUV] * ce;

			r_c=1;
			for (alpha = 0; alpha<=360; alpha=alpha + alpha_step){

				x_1 = x_c + R*sin(alpha*PI/180);
				y_1 = y_c + R*cos(alpha*PI/180);

				x_2 = x_c + R*sin((alpha + alpha_step)*PI/180);
				y_2 = y_c + R*cos((alpha + alpha_step)*PI/180);

				makePolygon(x_c/pixelWidth,y_c/pixelWidth,x_1/pixelWidth,y_1/pixelWidth,x_2/pixelWidth,y_2/pixelWidth);

				roiManager("Add");
				roiManager("Select", newArray(0,r_c));
				
				roiManager("AND");
				roiManager("Add");
				roiManager("Deselect");
				roiManager("Select", r_c);
				roiManager("Delete");

				r_c = r_c + 1;
				

		
			}
		

		run("Clear Results");


		nROIs = roiManager("Count")-1;

		chan = 0;
		//-----------------channel 1
		
		if (ch1 == true){
			selectWindow (title0);
			roiManager("Select", 0);
			
			run("Set Measurements...", "mean min integrated redirect=None decimal=3");
			run("Measure");
			ch1_average = getResult("Mean", 0);

			ch1_average_corr = ch1_average/background[chan] - 1;

			if (ch1_average_corr < 0 || ch1_average_corr != ch1_average_corr) ch1_average_corr = 0;

			run("Clear Results");

			for (i=0; i<nROIs; i++){
 
				roiManager("Select", i+1);
				
				roiManager("Rename", i*alpha_step);

				run("Set Measurements...", "mean min redirect=None decimal=3");
				run("Measure");

				Max_intensity_membrane [i] = getResult("Mean", i)/background[chan] - 1;

				if (Max_intensity_membrane [i]!=Max_intensity_membrane [i] || Max_intensity_membrane [i] < 0){
					Max_intensity_membrane [i] = 0;
				}
			
				Angle_array [i] = i*alpha_step;
			}


			run("Clear Results");
			chan ++;
		}
//-----------channel 2----------------------

		if (ch2 == true){
			selectWindow (title1);
			roiManager("Select", 0);
			
			run("Set Measurements...", "mean min integrated redirect=None decimal=3");
			run("Measure");
			ch2_average = getResult("Mean", 0);
			
			ch2_average_corr = ch2_average/background[chan] - 1;

			if (ch2_average_corr < 0 || ch2_average_corr != ch2_average_corr) ch2_average_corr = 0;
			
			run("Clear Results");

			for (i=0; i<nROIs; i++){
 
				roiManager("Select", i+1);
				roiManager("Rename", i*alpha_step);

				run("Set Measurements...", "mean min redirect=None decimal=3");
				run("Measure");

			

				Max_intensity_protein_1 [i] = getResult("Mean", i)/background[chan] - 1;
					if (Max_intensity_protein_1 [i]!=Max_intensity_protein_1 [i] || Max_intensity_protein_1 [i] < 0){
						Max_intensity_protein_1 [i] = 0;
					}
				
				Angle_array [i] = i*alpha_step;
			}


			run("Clear Results");
			chan ++;
		}
		
		//-----------channel 3----------------------
	
		if (ch3 == true){
			selectWindow (title2);
			roiManager("Select", 0);
			
			
			run("Set Measurements...", "mean min integrated redirect=None decimal=3");
			run("Measure");
			ch3_average = getResult("Mean", 0);
			

			ch3_average_corr = ch3_average/background[chan] - 1;

			if (ch3_average_corr < 0 || ch3_average_corr != ch3_average_corr) ch3_average_corr = 0;
			
			run("Clear Results");

			for (i=0; i<nROIs; i++){
 
				roiManager("Select", i+1);
				
				roiManager("Rename", i*alpha_step);

				run("Set Measurements...", "mean min redirect=None decimal=3");
				run("Measure");

			

				Max_intensity_protein_2 [i] = getResult("Mean", i)/background[chan] - 1;
					if (Max_intensity_protein_2 [i]!=Max_intensity_protein_2 [i] || Max_intensity_protein_2 [i] < 0){
						Max_intensity_protein_2 [i] = 0;
					}
				
				Angle_array [i] = i*alpha_step;
			}


			run("Clear Results");
			chan ++;
		}		
		
//-----------------------------------Separate profiles data output-----------------------------------------
//---------------------------------------------------------------------------------------------------------
		
		 xl=File.open(file_excell + "_GUV_"+number_of_GUV+1+".xls");
		
		
		print (xl, "PROTEIN VALUE: <protein>/background");
		print (xl, " ");
		print (xl, "If PROTEIN VALUE = NaN it goes to 0");
		print (xl, " ");

		row = "\t\t" + "angle";

		if (ch1 == true){
			row = row +"\t\t" + title0;
		}

		if (ch2 == true){
			row = row +"\t\t" + title1;
		}

		if (ch3 == true){
			row = row +"\t\t" + title2;
		}

		print(xl,row);
		
		

		for(i=0;i<nROIs;i++){

			row = "\t\t" + Angle_array[i];

		if (ch1 == true){
			row = row +"\t\t" + Max_intensity_membrane[i];
		}

		if (ch2 == true){
			row = row +"\t\t" + Max_intensity_protein_1[i];
		}

		if (ch3 == true){
			row = row +"\t\t" + Max_intensity_protein_2[i];
		}

		print(xl,row);

			
			
		}
		File.close(xl);	


//--------------------------------------------------------------------------------------------------------
//--------------------------------------Preview output----------------------------------------------------

if (preview_mode == "First file" || preview_mode == "Select file"){
	plot_name = "GUV " + number_of_GUV+1;

	Array.getStatistics(Max_intensity_membrane, min_1, max_1, mean_1, std_1);
	Array.getStatistics(Max_intensity_protein_1, min_2, max_2, mean_2, std_2);
	Array.getStatistics(Max_intensity_protein_2, min_3, max_3, mean_3, std_3);

	max_max = max_1;
	if (max_2 > max_max) max_max = max_2;
	if (max_3 > max_max) max_max = max_3;

	legend = "";
	
	
	Plot.create(plot_name, "angle, deg", "Intensity, a.u.");
	Plot.setLimits (0,360,-1,max_max+1);
	
	if (ch1 == true){
	Plot.setColor ("green");
	Plot.add ("Lines",Angle_array,Max_intensity_membrane);
	legend = legend + title0 + "\t";}

	if (ch2 == true){
	Plot.setColor ("red");
	Plot.add ("Lines",Angle_array,Max_intensity_protein_1);
	legend = legend + title1 + "\t";}

	if (ch3 == true){
	Plot.setColor ("blue");
	Plot.add ("Lines",Angle_array,Max_intensity_protein_2);
	legend = legend + title2 + "\t";}
	Plot.setColor ("green");

Plot.setLegend(legend, "top-right");
Plot.show();

	
	file_plot = dir1 + File.separator + list[number_of_file] +"_GUV" + number_of_GUV+1 + "_plot";
		

		saveAs("Tiff", file_plot);
		close();
}



		
//---------------------------------------------------------------------------------------------------------
		
//-----------------------------------Construct rows for summary files--------------------------------------
//---------------------------------------------------------------------------------------------------------
		chan = 0;

		if (ch1 == true){
			summary_ch1_unit  = list [number_of_file]+" GUV_" + number_of_GUV+1 + "\t\t"+ ch1_average + "\t\t"+ ch1_average_corr + "\t\t" + background[chan];
			chan++;
			summary_ch1 = Array.concat(summary_ch1,summary_ch1_unit);
		}
		
	    if (ch2 == true){
			summary_ch2_unit  = list [number_of_file]+" GUV_" + number_of_GUV+1 + "\t\t"+ ch2_average + "\t\t"+ ch2_average_corr + "\t\t" + background[chan];
			chan++;
			summary_ch2 = Array.concat(summary_ch2,summary_ch2_unit);
		}
	    if (ch3 == true){
			summary_ch3_unit  = list [number_of_file]+" GUV_" + number_of_GUV+1 + "\t\t"+ ch3_average + "\t\t"+ ch3_average_corr + "\t\t" + background[chan];
			chan++;
			summary_ch3 = Array.concat(summary_ch3,summary_ch3_unit);
		}
		}


run("Close All");
		
		
		
		if (preview_mode == "First file" || preview_mode == "Select file") {
		//open(file_marked);
		Prename=dir1 + File.separator + list[number_of_file] +"_GUV";
		number_of_file = list.length-1;
		
		
		}
	
	
	
	}
	
	}

	ScreenClean();


if (preview_mode == "First file" || preview_mode == "Select file") {
		open(file_marked);
		for (number_of_GUV = 0; number_of_GUV < GUVs_total; number_of_GUV ++){
			file_plot = Prename + number_of_GUV+1 + "_plot.tif";
			open(file_plot);
		}
		
		}




//-----------------------------------Summary files output--------------------------------------------------
//---------------------------------------------------------------------------------------------------------
if (ch1 ==true){
xl1=File.open(dir1+ File.separator+"Summary_" + title0 + ".xls");

print(xl1,"name"+"\t\t"+  "<"+title0+">" +"\t\t" + "<"+title0+">/background" +"\t\t" + "background");

for (number_of_file = 0; number_of_file<summary_ch1.length; number_of_file++){
	
		print(xl1,summary_ch1[number_of_file]);	
	
}

File.close(xl1);
}

if (ch2 ==true){
xl1=File.open(dir1+ File.separator+"Summary_" + title1 + ".xls");

print(xl1,"name"+"\t\t"+  "<"+title1+">" +"\t\t" + "<"+title1+">/background" +"\t\t" + "background");

for (number_of_file = 0; number_of_file<summary_ch2.length; number_of_file++){
	
		print(xl1,summary_ch2[number_of_file]);	
	
}

File.close(xl1);
}

if (ch3 ==true){
xl1=File.open(dir1+File.separator+"Summary_" + title2 + ".xls");

print(xl1,"name"+"\t\t"+  "<"+title2+">" +"\t\t" + "<"+title2+">/background" +"\t\t" + "background");

for (number_of_file = 0; number_of_file<summary_ch3.length; number_of_file++){
	
		print(xl1,summary_ch3[number_of_file]);	
	
}

File.close(xl1);
}
//---------------------------------------------------------------------------------------------------------


}








// ==========================================================================================================
// ***********************************            Functions         *****************************************
// ==========================================================================================================

//-----------------------------------Delete array elements (array, array of indexes to delete)-------------
//---------------------------------------------------------------------------------------------------------
function delete_array_elements (arr , index_arr){
	Array.sort(index_arr);
	
	subt = 0;
	for (ind = 0; ind<index_arr.length; ind++){
		arr = delete_array_element (arr, index_arr[ind]+1 - subt);
		subt++;
	}
		return arr;
}
//---------------------------------------------------------------------------------------------------------

//-----------------------------------Delete array element (array,index to delete)--------------------------
//---------------------------------------------------------------------------------------------------------
function delete_array_element (arr , index){
	a = Array.slice(arr,index);
	
	b = Array.slice(arr,0,index-1);
	c = Array.concat (b,a);
	return c;
}
//---------------------------------------------------------------------------------------------------------

//---------------Get coordinates of all points in an active ROI and save them in an array------------------
//---------------------------------------------------------------------------------------------------------
function get_all_points_in_a_ROI (){
	Stack.getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);

	points = newArray;
	

	
	for (x = 0; x<=width; x++){
		for (y = 0; y<=height; y++){
			if (Roi.contains(x, y)==true){
				points = Array.concat (points,x);
				points = Array.concat (points,y);
			}
		}
	}

	return points;
}
//---------------------------------------------------------------------------------------------------------

//---------------Get list of x coordinates of all points in an active ROI----------------------------------
//---------------------------------------------------------------------------------------------------------
function get_x(arr){
	arr1 = newArray;

	for (i=0; i<arr.length; i+=2)
	arr1 = Array.concat (arr1,arr[i]);

	return arr1;
}
//---------------------------------------------------------------------------------------------------------

//---------------Get list of y coordinates of all points in an active ROI----------------------------------
//---------------------------------------------------------------------------------------------------------
function get_y(arr){
	arr1 = newArray;

	for (i=1; i<arr.length; i+=2)
	arr1 = Array.concat (arr1,arr[i]);

	return arr1;
}
//---------------------------------------------------------------------------------------------------------

//---------------Clean screen from all opened windows------------------------------------------------------
//---------------------------------------------------------------------------------------------------------
function ScreenClean()
      {	
	while (nImages>0) close();

          WinOp=getList("window.titles");
	for(i=0; i<WinOp.length; i++)
	  {selectWindow(WinOp[i]);run ("Close");}

	  fenetres=newArray("B&C","Channels","Threshold");
	for(i=0;i!=fenetres.length;i++)
	   if (isOpen(fenetres[i]))
	    {selectWindow(fenetres[i]);run("Close");}
       }
//---------------------------------------------------------------------------------------------------------

