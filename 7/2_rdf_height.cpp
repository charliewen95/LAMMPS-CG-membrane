#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
#include<vector>
#include<cstdlib>
using namespace std;

struct Bead {
    int id;
    int type;
    int mol;
    double x;
    double y;
    double z;
};
	
double dist_form_3d(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double dist = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));
    return dist;
}

double dist_form_2d(double x1, double y1, double x2, double y2)
{
    double dist = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
    return dist;
}

int main() {
    string fileName1;
    ifstream file1;
    cout << "Which data file to process?";
    cin >> fileName1;
    file1.open(fileName1);
    ofstream output1;
    output1.open("height_rdf.out");

    string line;
    Bead bd;
    bool bb=true;
    vector<double> b3tob4;
    double a,b,c,d,e,f,rdist;
    double centerprot[3]={-8.331560134887695,6.462148666381836,3.196976661682129};
    
    output1 << "mol height distance" << endl;
    
    if (file1.is_open()) {
        while(getline(file1,line)) {
            istringstream iss(line);
            if (iss >> a >> b >> c >> d >> e >> f) {
                bd.id=a;
                bd.type=b;
                bd.mol=c;
                bd.x=d;
                bd.y=e;
                bd.z=f;
                rdist=dist_form_2d(bd.x,bd.y,centerprot[0],centerprot[1]);
                output1 << bd.mol << " " << bd.z << " " << rdist << endl;
            } else {
                cout << "wrong format." << endl;
            }
        }
    }
}

