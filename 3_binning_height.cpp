#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
#include<vector>
#include<cstdlib>
#include<limits>
using namespace std;
double min_max(string fileName) {
    ifstream inputFile1;
    inputFile1.open(fileName);
    if (!inputFile1.is_open()) {
        cerr << "Error opening file " << fileName << endl;
        return 1;
    }

    double a{}, b{}, c{}, d{};
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();

    if (inputFile1.is_open()) {
        string line;
        while(getline(inputFile1, line)) {
            istringstream iss(line);
            if (iss >> a >> b >> c >> d) {
                if (d < min) {
                    min = d;
                }
                if (d > max) {
                    max = d;
                }
            }
        }
    }
    cout << "min: " << min << '\n';
    cout << "max: " << max << '\n';
    return min,max;
}

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
    outputFile1.open("height_binned.out");
    if (!outputFile1.is_open()) {
        cerr << "Error opening output file." << endl;
        return 2;
    }

    string line;
    int lineCount = 0;

    int numOfBins=190;
    int numOfDatToBin=3;
    //double low,high=min_max(fileName);
    double low=0;
    //double binWidth=(836.853-0.484321)/100;
    double binWidth=5;

    // 2D vector -->  vector<vecotr<double>> binning(rows, vector<double>(columns))
    // 2D vector -->  vector<vecotr<double>> binning(number_of_bins, vector<double>(data_to_be_binned))
    vector<vector<double>> binning(numOfBins, vector<double>(numOfDatToBin));

    while (getline(inputFile1, line)) {
        double a, b, c;
        istringstream iss(line);
        if (iss >> a >> b >> c) {
            double binNum = floor(c/binWidth);
            binning[binNum][0]++;
            binning[binNum][1]+=b;
        } else {
            cerr << "Error parsing line " << lineCount + 1 << endl;
            continue;
        } 
    }
    outputFile1 << "height rdist" << endl;
    for (int i=0;i<numOfBins;i++){
        double label=i*binWidth+low+(binWidth/2);
        outputFile1 << (binning[i][1]/binning[i][0]) << " " << label << endl; 
    }
    ////////// END //////////
    cout << "DONE" << endl;
    inputFile1.close();
    outputFile1.close();
}
