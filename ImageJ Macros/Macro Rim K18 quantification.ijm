// open stack
// duplicate single slice
// close stack

//create masks of K18 and DAPI channel
rename("disc");
run("Split Channels");
selectWindow("C2-disc");
run("Duplicate...", "title=K18 duplicate");
setAutoThreshold("Default dark");
//run("Threshold...");
setThreshold(10, 255);
run("Create Mask");
rename("mask_K18");
run("Fill Holes");
selectWindow("K18")
close()
selectWindow("C1-disc");
run("Duplicate...", "title=dapi duplicate");
setAutoThreshold("Default dark");
//run("Threshold...");
setThreshold(10, 255);
run("Create Mask");
rename("mask_dapi");
run("Fill Holes");
selectWindow("dapi")
close()

//create a combined mask and an eroded combined mask
imageCalculator("Add create", "mask_dapi","mask_K18");
selectWindow("Result of mask_dapi");
rename("mask_combined");
run("Duplicate...", "title=mask_combined_eroded");
run("Options...", "iterations=10 count=1 black do=Nothing");
run("Erode");

// measure integrated intensity of K18 
selectWindow("mask_combined_eroded");
run("Divide...", "value=255");
imageCalculator("Multiply create", "C2-disc","mask_combined_eroded");
selectWindow("Result of C2-disc");
rename("K18_eroded");
run("Set Measurements...", "integrated redirect=None decimal=6");
run("Measure");
selectWindow("mask_combined");
run("Divide...", "value=255");
imageCalculator("Multiply create", "C2-disc","mask_combined");
selectWindow("Result of C2-disc");
rename("K18_whole");
run("Measure");

// measure integrated intensity of mask as a measure of area (each pixel = 1 intensity) (for normalization)
selectWindow("mask_combined_eroded");
run("Measure");
selectWindow("mask_combined");
run("Measure");

//measure DAPI (integrated intensity for normalization)
imageCalculator("Multiply create", "C1-disc","mask_combined_eroded");
selectWindow("Result of C1-disc");
rename("dapi_eroded");
run("Measure");
imageCalculator("Multiply create", "C1-disc","mask_combined");
selectWindow("Result of C1-disc");
rename("dapi_whole");
run("Measure");

//close windows
selectWindow("mask_K18");
close();
selectWindow("mask_dapi");
close();
selectWindow("mask_combined");
close();
selectWindow("mask_combined_eroded");
close();
selectWindow("C1-disc");
close();
selectWindow("C2-disc");
close();

//order of measurements
//intensity of H3K18ac of eroded mask/disc
//intensity of H3K18ac of whole mask/disc
//"area" of eroded mask/disc
//"area" of whole mask/disc
//intensity of DAPI of eroded mask/disc
//intensity of DAPI of whole mask/disc
//=> use RawIntDen!!