// Measures nuclear vs total red signal by using DAPI  channel  (blue) as a mask
// Need to first draw a region of interest in the image before running the macro
// The first row is total signal,  the 2nd row is nuclear

run("Set Measurements...", "integrated redirect=None decimal=3");

if(roiManager("count")>0){
	roiManager("deselect");
	roiManager("delete");
}
roiManager("Add");

rename("disc");
run("Split Channels");

selectWindow("disc (green)");
close();

selectWindow("disc (blue)");
setAutoThreshold("Default dark");
run("Convert to Mask");
run("Divide...", "value=255");

selectWindow("disc (red)");
run("Duplicate...", " ");
imageCalculator("Multiply", "disc (red)-1","disc (blue)");

selectWindow("disc (red)");
roiManager("Select", 0);
run("Measure");

selectWindow("disc (red)-1");
roiManager("Select", 0);
run("Measure");
