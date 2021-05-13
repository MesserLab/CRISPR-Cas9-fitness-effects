//Anna Maria Langmueller for the Messer Lab
//February 2020 


//ImageJ macro to classifiy RR, RG, and GG flies on Screen 
//Input: .JPG images (w, r, g : in that order)
//Program will detect flies, zoom into every single fly, and parse through the three pictures
//to change speed of picture parsing, change the value on line 232 (wait(330)) 
//User has to classify manually 
//make sure you ticked "Add to Manager" in measurements!
//requires Template_Matching Plugin: https://sites.google.com/site/qingzongtseng/template-matching-ij-plugin#downloads 

//////////////////
//FUNCTIONS
//////////////////

//function clean_up: 
//close all previous images and remove entries in  ROI Manager

function clean_up() {
	run("Clear Results");
	roiManager("reset");
	while (nImages>0){
		selectImage(nImages);
		close();
	}
}

//function filter_ambient:
//applies filtering on white light picture to obtain outlines of flies

function filter_ambient(my_ambient_id){
	selectImage(my_ambient_id);
	run("Median...", "radius=3");
	run("Gaussian Blur...", "sigma=3");
	run("Auto Threshold", "method=MaxEntropy ignore_black white");
	run("Fill Holes");
	run("Open");
	run("Analyze Particles...", "size=750-Infinity show=Overlay display exclude include add");
}

//function enlarge_flies:
//enlarge convex hull of each fly by 15 pixels and store bigger outline

function enlarge_flies(){
	numROIS=roiManager("count");
	print(numROIS);
	for(i=0; i<numROIS; i++) {
		roiManager("select",i);
		roiManager("rename",i+1);
		run("Convex Hull");
		run("Enlarge...","enlarge=20");
		roiManager("update");
		roiManager("deselect");
		}
run("Select None");
numROIS=roiManager("count");
print(numROIS);
}

//function filter_fluorescent:
//applies median filter of 2 pixels to fluorescent picture

function filter_fluorescent(myID) {
	selectImage(myID);
	run("Median...", "radius=2");
}

//////////////////
//MACRO STEPS
//////////////////

//////////////////
//Preparation
//////////////////

//clean up ImageJ
clean_up();

//user input
//ask for the directory
ambient=File.openDialog("Select ambient light picture"); 
red=File.openDialog("Select fluorescent red picture");
green=File.openDialog("Select fluorescent green picture");

//open + store ID of raw pictures
open(ambient);
path_ambient=getDirectory("image");
print(path_ambient);
ambientID=getImageID();
open(red);
redID=getImageID();
path_red=getDirectory("image");
print(path_red);
open(green);
greenID=getImageID();
path_green=getDirectory("image");
print(path_green);

//////////////////
//Modifications
//////////////////

//modification of raw pictures
//white: splite channels, keep red and green 
selectImage(ambientID);
run("Split Channels");
ambient_blueID=getImageID();
ambient_greenID=ambient_blueID+1;
ambient_redID=ambient_blueID+2;
selectImage(ambient_blueID);
close();
//red: splict channels, keep only red 
selectImage(redID);
run("Split Channels");
red_blueID=getImageID();
red_greenID=red_blueID+1;
red_redID=red_blueID+2;
selectImage(red_blueID);
close();
selectImage(red_greenID);
close();
selectImage(red_redID);
filter_fluorescent(red_redID);
//green: split channels, keep only green 
selectImage(greenID);
run("Split Channels");
green_blueID=getImageID();
green_greenID=green_blueID+1;
green_redID=green_blueID+2;
selectImage(green_blueID);
close();
selectImage(green_redID);
close();
selectImage(green_greenID);
filter_fluorescent(green_greenID);

////////////////////////
//Stack-Crop-Alignment
////////////////////////

//remaining images: ambient_greenID, ambient_redID, red_redID, green_greenID
//stack images: 
run("Images to Stack", "name=alignment title=[] use");

//Dialog with user - do you need to crop?
Dialog.create("Cropping of Image");
Dialog.addChoice("Do you want to crop?",newArray("yes","no"),"no");
Dialog.show();
cropping=Dialog.getChoice();
if(cropping == "yes") { 
	waitForUser("Select Region");
	run("Crop");
	run("Select None");
}

//Dialog with user - do you want to align?
Dialog.create("Alignment Stack");
Dialog.addChoice("Do you want to align?",newArray("yes","no"),"yes");
Dialog.show();
alignment=Dialog.getChoice();
if(alignment== "yes") {
	//run slice alignment (needs user input)
	run("Align slices in stack...");
}

//split aligned stack back to pictures 
run("Stack to Images");
//new IDs neccessary (reverse order)
green_greenID=getImageID();
red_redID=green_greenID+1;
ambient_greenID=green_greenID+2;
ambient_redID=green_greenID+3;


////////////////////////
//Find Flies
////////////////////////

//substract ambient_green from ambient_red 
imageCalculator("subtract", ambient_redID, ambient_greenID);
ambient_modifiedID=ambient_redID;
selectImage(ambient_modifiedID);

//enhancement on ambient light picture 
Dialog.create("Maximum filter");
Dialog.addChoice("Run Renyi Entropy instead of Max Entropy?", newArray("yes","no"),"no");
Dialog.show();
myMax=Dialog.getChoice();
if(myMax=="yes") {
	selectImage(ambient_modifiedID);
	run("Median...", "radius=3");
	run("Gaussian Blur...", "sigma=3");
	run("Auto Threshold", "method=RenyiEntropy ignore_black white");
	run("Fill Holes");
	run("Open");
	run("Analyze Particles...", "size=750-Infinity show=Overlay display exclude include add");
}

if(myMax=="no") {
//determine outlines of the flies 
filter_ambient(ambient_modifiedID);
}
//enlrage the outlines of the flies

enlarge_flies();
//close binary image
selectImage(ambient_modifiedID);
close();

////////////////////////
//Zoom to single Flies
////////////////////////

run("Images to Stack", "name=Stack");
run("From ROI Manager");
numROIS=roiManager("count");
print(numROIS);


Stack.setSlice(1);
//select the Stack 
selectWindow("Stack");
//for each ROI

for(i=0; i<numROIS; i++) {
	//zoom into a certain selection 
	roiManager("select",i);
	run("To Selection");
	for (j = 0; j < nSlices; j++) {
		Stack.setSlice(j + 1);
		wait(330);
	}
}

run("Select None");
numROIS=roiManager("count");
print(numROIS);
