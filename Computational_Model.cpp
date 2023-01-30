#include <sys/stat.h>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <ctime>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <random>
#include <algorithm>

//THIS CODE REQUIRES C++11

class Exciton{
public:
    std::string name;
    int X_pos;
    int Y_pos;
    bool active;
    int decayX;
    int decayY;
    std::string fate;
    int decayTime;

    Exciton(int X_pos, int Y_pos, bool active, std::string name){
        this->X_pos = X_pos;
        this->Y_pos = Y_pos;
        this->active = active;
        this->name = name;
        decayTime=0;
    }
};

//GLOBAL
using std::cout;
using std::endl;
using std::string;
using std::vector;
int spawn_index = 0;

//Define the active exciton list. This is what will keep track of the "living" exciton population
vector<Exciton> active_Excitons;

//Set spawn seet for exciton spawning (we only need 1 seed for the spawn gaussian each time we run the program so I put it in global)
std::random_device rd;
std::mt19937 rng(rd());
std::uniform_real_distribution<double> uni(0,1);

//Nx and Ny determine the size of the grid
const int Nx = 100;
const int Ny = 100;

//initialize functions
double bivariate_gaussian(double x, double y);
void spawn_exciton(double probability, int x, int y);
bool decay_chance(int cycle_index, Exciton exciton);
int randomNum1_4(void);
Exciton move_exciton(int directionNum, Exciton exciton);


//define the cell struct to store x and y points which we will access like ordered pairs
struct cell{
    int x;
    int y;
};

int main(){
    //CREATE FOLDER FOR THE SIMULATION OUTPUT
    // Get current date and time
    std::time_t t = std::time(nullptr);
    char date_time[100];
    std::strftime(date_time, sizeof(date_time), "%Y-%m-%d %H:%M:%S", std::localtime(&t));

    // Create directory path with unique name
    string custom_path = "/Users/david/Desktop/C++/Exciton Research/Heat Map R/" + string(date_time) + string(" Data");

    // Check if the directory already exists or not
    struct stat st;
    if (stat(custom_path.c_str(), &st) == 0) {
        cout << "Error: Directory already exists \n";
    } else {
        // Create directory
        int status = mkdir(custom_path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        if (status == 0) {
            cout << "Directory created successfully at " << custom_path << "\n";
        } else {
            cout << "Error creating directory: " << strerror(errno) << "\n";
        }
    }

    //initialize the cell grid
    cell grid[Nx][Ny];

    //Cycle through every point on the grid (INITIALIZATION LOOP) this loop runs only once
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            //Grid points take the form (i,j) and we're assigning numerical values to the ordered pairs stored in the matrix
            grid[i][j].x = i;
            grid[i][j].y = j;

            //Default spawn probability is 0
            double Gauss_height = 0;
            //grab height of gaussian (used to get spawn probability)
            Gauss_height = bivariate_gaussian(grid[i][j].x, grid[i][j].y);
            
            //determine whether or not an exciton spawns based on the statistical distribution (from bivariate gaussian)
            spawn_exciton(Gauss_height, grid[i][j].x, grid[i][j].y);
        }   
    }
// Create the file name and path
string file_path = custom_path + "/Excitons_0.csv";

std::ofstream excitonFile;
excitonFile.open(file_path);

//Write headers for columns
excitonFile << "X_pos, Y_pos, \n";

//Cycle through each exciton in the active exciton list
for (int i=0; i<active_Excitons.size(); i++){
    //Write values for each exciton into the appropriate column
    excitonFile << active_Excitons.at(i).X_pos << ", " << active_Excitons.at(i).Y_pos << "\n";
}
excitonFile.close();
int index1 = 1;

    //MAIN LOOP
    int j=0;
    while (j<4000){
        //Chance to decay randomly

        //check if other excitons are nearby and run annihilation function if so

        //Get the total # of living excitons
        int exciton_population_size = active_Excitons.size();

        //Cycle through all living excitons
        for (int i = 0; i < exciton_population_size; i++) {
            //Check to see if the exciton has decayed, if so, don't run the time evolution code block
            if (active_Excitons.at(i).active){
                //Generate a random number between 1 and 4
                int randomDirection;
                randomDirection = randomNum1_4();

                //Move the exciton in a totally random direction
                active_Excitons.at(i) = move_exciton(randomDirection, active_Excitons.at(i));

                //Add 1 to the number of cycles the exciton has been alive
                active_Excitons.at(i).decayTime++;
            }
        }
        //Delete dead excitons from the system
        int index = 0;
        for (int i = 0; i < active_Excitons.size(); i++) {
            if (active_Excitons[i].active == true) {
                std::swap(active_Excitons[index], active_Excitons[i]);
                index++;
            }
        }
        active_Excitons.erase(active_Excitons.begin() + index, active_Excitons.end());
        if (j%25 == 0){
        //DATA STORAGE
        // Create a new csv file for the exciton locations
        string file_path = custom_path + "/Excitons_";
        //Add a 1 to the end of the data name for each output (impractical temporary solution to Rstudio bug)
        for (int i=0; i<index1; i++){
            file_path += std::to_string(1);
        }
        file_path += std::to_string(1) + ".csv"; 
        std::ofstream excitonFile;
        excitonFile.open(file_path);

        // Write headers for columns
        excitonFile << "X_pos, Y_pos, \n";

        // Cycle through each exciton in the active exciton list
        for (int i=0; i<active_Excitons.size(); i++){
            // Write values for each exciton into the appropriate column
            excitonFile << active_Excitons.at(i).X_pos << ", " << active_Excitons.at(i).Y_pos << "\n";
        }

        excitonFile.close();
        index1++;
        }
        //increase while-loop index
        j++;
    }

}
//This function generates our gaussian distribution stamp
double bivariate_gaussian(double x, double y) {
    //Bivariate gaussian parameters
    const double A = 1.0; // amplitude (if A =/= 1 then statistics will break)
    const double x0 = (Nx-1)/2.0; // center x -- default to halfway along x grid
    const double y0 = (Ny-1)/2.0; // center y -- default to halfway along x grid (this centers the laser impulse at the center of the grid)
    const double sigma_x = Nx / 4.0; // standard deviation along x-axis
    const double sigma_y = Ny / 4.0; // standard deviation along y-axis
    double exponent = -((x - x0) * (x - x0) / (2 * sigma_x * sigma_x) + (y - y0) * (y - y0) / (2 * sigma_y * sigma_y));

    //return the height of the bivariate gaussian at the point (x,y)
    return A * std::exp(exponent);
}

//This function spawns excitons based onto the grid based on the probability from the gaussian
void spawn_exciton(double probability, int x, int y){
    //If exciton amplitude of gaussian is not 1 return this error  message
    if (probability < 0 || probability > 1) {
        cout << "Invalid probability input, probability must be between 0 and 1 \n";
        return;
    }

    //Randomly generate a number between 0 and 1
    double randomNumber = uni(rng);

    if (randomNumber < probability) {
        //Generate exciton label
        std::string exciton_name = "exciton_";
        std::string exciton_index = std::to_string(spawn_index);
        exciton_name += exciton_index;
        spawn_index++;

        //Instantiate unique exciton object using the label we just generated
        Exciton exciton(x,y,true,exciton_name);

        //Add the exciton to the list of active (not decayed) excitons
        active_Excitons.push_back(exciton);
    }
}

//Generate a random integer between 1 and 4
int randomNum1_4(void){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(1, 4);
    return dist(mt);
}

Exciton move_exciton(int directionNum, Exciton exciton){
    //returns an error if the random number used to pick a direction is not between 1 and 4 or not an integer
    //Commented out for efficiency
    // if (directionNum != 1 && directionNum != 2 && directionNum != 3 && directionNum != 4){
    //     cout<<"Error: Did not recieve a direction \n";
    // }
    if (directionNum==1){
        exciton.X_pos += 1;
        //Check to make sure exciton is inside the grid
        if (exciton.X_pos < Nx && exciton.X_pos > 0 && exciton.Y_pos < Ny && exciton.Y_pos > 0){
            return exciton;
        }
        //If exciton is not in grid, delete it
        else {
            exciton.active = false;
            exciton.fate = "Out of bounds, deletion";
            return exciton;
        }

    }
    else if (directionNum==2){
        exciton.Y_pos+=1;
        //Check to make sure exciton is inside the grid
        if (exciton.X_pos < Nx && exciton.X_pos > 0 && exciton.Y_pos < Ny && exciton.Y_pos > 0){
            return exciton;
        }
        //If exciton is not in grid, delete it
        else {
            exciton.active = false;
            exciton.fate = "Out of bounds, deletion";
            return exciton;
        }
    }
    else if (directionNum==3){
        exciton.X_pos-=1;
        //Check to make sure exciton is inside the grid
        if (exciton.X_pos < Nx && exciton.X_pos > 0 && exciton.Y_pos < Ny && exciton.Y_pos > 0){
            return exciton;
        }
        //If exciton is not in grid, delete it
        else {
            exciton.active = false;
            exciton.fate = "Out of bounds, deletion";
            return exciton;
        }
    }
    else if (directionNum==4){
        exciton.Y_pos-=1;
        //Check to make sure exciton is inside the grid
        if (exciton.X_pos < Nx && exciton.X_pos > 0 && exciton.Y_pos < Ny && exciton.Y_pos > 0){
            return exciton;
        }
        //If exciton is not in grid, delete it
        else {
            exciton.active = false;
            exciton.fate = "Out of bounds, deletion";
            return exciton;
        }
    }
    cout<<"move_exciton function failed to return after movemenet, check for bugs.\n";
    return exciton;
}


bool decay_chance(int cycle_index, Exciton exciton){
    //randomly generate a number between 1 and 100 (or whatever the spontaneous decay rate is going to be)
    int randDecayNum;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(1, 100);
    randDecayNum = dist(mt);

    if (randDecayNum == 1){
        exciton.active = false;
        exciton.decayX = exciton.X_pos;
        exciton.decayY = exciton.Y_pos;
        exciton.fate = "QE before localization";
        
        return true;
    }

    else{
        return false;
    }
}

 //TEST/DEBUG FUNCTIONS

    //PRINT OUT THE POSITION, NAME, AND ACTIVE STATE OF ALL EXCITONS IN THE active_Excitons VECTOR
    // for (int i=0; i<active_Excitons.size(); i++){
    //     cout<<"obj name is "<<active_Excitons.at(i).name<<"\n";
    //     cout<<"x_pos= "<< active_Excitons.at(i).X_pos<<"\n";
    //     cout<<"y_pos= "<<active_Excitons.at(i).Y_pos<<"\n";
    //     cout<<active_Excitons.at(i).active<<"\n";
    // }

    // PRINT OUT THE GRID IN ORDERED PAIRS
    // for (int i = 0; i < Ny; i++) {
    //     for (int j = 0; j < Nx; j++) {
    //         std::cout << "(" << grid[j][i].x << "," << grid[j][i].y << ") ";
    //     }
    //     std::cout << std::endl;
    // }

    //PRINT BIVARIATE GAUSSIAN OUTPUTS
    // double test;
    // test = bivariate_gaussian(grid[5][5].x, grid[5][5].y);
    // cout<<test<<endl;

    //-------------

    //TEST RANDOM NUMBER GENERATOR
    // for (int i = 0; i < 10; i++) {
    // int random_num = randomNum1_4();
    // cout << random_num << endl;
    //}

    //TEST MAIN LOOP (paste this below the main loop to follow the activity of an exciton before and after the main loop)
    // cout<<"obj name is "<<active_Excitons.at(7).name<<"\n";
    // cout<<"x_pos= "<< active_Excitons.at(7).X_pos<<"\n";
    // cout<<"y_pos= "<<active_Excitons.at(7).Y_pos<<"\n";
    // cout<<"decay time= "<<active_Excitons.at(7).decayTime<<"\n";

        //TEST EXCITON MOVEMENT
    // int randomDirection = randomNum1_4();

    // Exciton myExciton(50,50,true,"TestExciton");
    // active_Excitons.push_back(myExciton);

    // for (int i=0; i<active_Excitons.size(); i++){
    //     cout<<"obj name is "<<active_Excitons.at(i).name<<"\n";
    //     cout<<"x_pos= "<< active_Excitons.at(i).X_pos<<"\n";
    //     cout<<"y_pos= "<<active_Excitons.at(i).Y_pos<<"\n";
    //     cout<<active_Excitons.at(i).active<<"\n";
    // }

    // active_Excitons.at(0) = move_exciton(randomDirection, active_Excitons.at(0));

    // for (int i=0; i<active_Excitons.size(); i++){
    //     cout<<"obj name is "<<active_Excitons.at(i).name<<"\n";
    //     cout<<"x_pos= "<< active_Excitons.at(i).X_pos<<"\n";
    //     cout<<"y_pos= "<<active_Excitons.at(i).Y_pos<<"\n";
    //     cout<<active_Excitons.at(i).active<<"\n";
    // }

