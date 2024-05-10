#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
#include<vector>
#include<cstdlib>
#include<limits>
using namespace std;
int main() {
    string fileName;
    cout << "Enter the name of the data file: ";
    cin >> fileName;

    ifstream inputFile1;
    inputFile1.open(fileName);
    if (!inputFile1.is_open()) {
        cerr << "Error opening file " << fileName << endl;
        return 1;
    }
    
    // std::ios::out tells the program to open the file in output mode, which means that we can write to the file.
    // std::ios::app tells the program to open the file in "append" mode, which means that any data written to the file will be appended to the end of the file. If the file does not exist, it will be created.
    // ofstream outputFile1("binned.out", std::ios::out | std::ios::app);
    ofstream outputFile1;
    outputFile1.open("extract.out");
    if (!outputFile1.is_open()) {
        cerr << "Error opening output file." << endl;
        return 2;
    }

    string line;
    int lineCount = 0;

    while (getline(inputFile1, line)) {
        //if (line.find("ITEM: TIMESTEP")==string::npos) {
        //   cout<<line<<endl;
        //}
        double id,type,mol,x,y,z;
        istringstream iss(line);
        if (iss>>id>>type>>mol>>x>>y>>z) {
            if (type==2&&mol!=99999){//////////INSERT SEARCH COMMAND HERE
                outputFile1<<id<<" "<<type<<" "<<mol<<" "<<x<<" "<<y<<" "<<z<<endl;
            }
        } else {
            cerr << "Error parsing line " << lineCount + 1 << endl;
            continue;
        } 
    }
//    outputFile1 << "Column Title" << endl;
//    for (int i=0;i<numOfBins;i++){
//        outputFile1 << "DATA1" << " " << "DATA2" << endl; 
//    }
    ////////// END //////////
    cout << "DONE" << endl;
    inputFile1.close();
    outputFile1.close();
}
