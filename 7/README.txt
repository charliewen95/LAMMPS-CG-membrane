1_SEARCH_EXTRACT
--> extracts coordinate of all intermediate beads (Type 2) and not resid 99999 (PIEZO)

2_rdf_height
--> extracts "resid coordinate_z distance_from_center"
Note::: if there is a change of center of protein open 2_rdf_height.cpp and change value of
double centerprot[3]={-8.331560134887695,6.462148666381836,3.196976661682129};

3_binning_height
--> bins the data change in file if necessary
int numOfBins=190;
double binWidth=5;

before run

make 1_SEARCH_EXTRACT
make 2_rdf_height
make 3_binning_height