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
//Time not synced properly for real modelling this is for testing purposes

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
    double expectation_x_term1; //for <x^2>
    double expectation_x_term2; //for <x>^2
    double expectation_y_term1; //for <y^2>
    double expectation_y_term2; //for <y>^2
    

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

//Global strain attraction parameter
const double K_sA = 10;

//Turn on/off whether or not you want to count excitons in region1 (false=off)
bool region1 = true;

//Turn on/off dimensionless units (false = off)
bool dimensionless = true;

//Define the active exciton list. This is what will keep track of the "living" exciton population
vector<Exciton> active_Excitons;

//Set spawn seet for exciton spawning (we only need 1 seed for the spawn gaussian each time we run the program so I put it in global)
std::random_device rd;
std::mt19937 rng(rd());
std::uniform_real_distribution<double> uni(0,1);

//Nx and Ny determine the size of the grid
const int Nx = 100;
const int Ny = 100;

//define the cell struct to store x and y points which we will access like ordered pairs
struct cell{
    int x;
    int y;
    //Is this tile strained (set strained tiles manually)
    bool strained;
    double P_left;
    double P_right;
    double P_up;
    double P_down;
    //Parameter for width of gaussian potential well (manually set the width of the well in the initialization loop)
    double sigma_strained;
    //Parameter for depth(amplitude) of gaussian potential well (manually set in the initialization loop)
    double amplitude;

    //Total force and total force componenets for probability calculations
    double F_net;
    double F_net_x;
    double F_net_y;
};

//Define struct to store variances along x and y
struct variance{
    double x_variance_struct;
    double y_variance_struct;
};



//initialize functions
double bivariate_gaussian(double x, double y);
void spawn_exciton(double probability, int x, int y);
Exciton decay_chance(Exciton exciton);
int randomNum1_4(void);
Exciton move_exciton(int directionNum, Exciton exciton);
int weightedRandomNum1_4(double P_right, double P_up, double P_left, double P_down);
Exciton get_radii(Exciton exciton);
variance calculate_variance(variance variance_total);
double exciton_density_region_analysis1(void);


int main(){
    
    //CREATE FOLDERS FOR THE SIMULATION OUTPUT
    // Get current date and time
    std::time_t t = std::time(nullptr);
    char date_time[100];
    std::strftime(date_time, sizeof(date_time), "%Y-%m-%d %H:%M:%S", std::localtime(&t));

    // Create directory path with unique name
    string positons_path = "/Users/david/Desktop/C++/Exciton Research/Heat Map R/" + string(date_time) + string(" Positons");
    string variance_path = "/Users/david/Desktop/C++/Exciton Research/Heat Map R/" + string(date_time) + string(" Variance");
    string regional_density_path1 = "/Users/david/Desktop/C++/Exciton Research/Heat Map R/" + string(date_time) + string(" eDensity1");

    // Check if the directory already exists or not
    struct stat st;
    if (stat(positons_path.c_str(), &st) == 0) {
        cout << "Error: Positions directory already exists \n";
    } else {
        // Create position data directory
        int status = mkdir(positons_path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        if (status == 0) {
            cout << "Directory created successfully at " << positons_path << "\n";
        } else {
            cout << "Error creating directory: " << strerror(errno) << "\n";
        }
    }
    //Check if directory exists or not
    if (stat(variance_path.c_str(), &st) == 0) {
        cout << "Error: Variance directory already exists \n";
    } else {
        // Create variance data directory
        int status = mkdir(variance_path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        if (status == 0) {
            cout << "Directory created successfully at " << variance_path << "\n";
        } else {
            cout << "Error creating directory: " << strerror(errno) << "\n";
        }
    }


    //Check if directory exists or not
    if (stat(regional_density_path1.c_str(), &st) == 0) {
        cout << "Error: Variance directory already exists \n";
    } else {
        // Create variance data directory
        int status = mkdir(regional_density_path1.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        if (status == 0) {
            cout << "Directory created successfully at " << regional_density_path1 << "\n";
        } else {
            cout << "Error creating directory: " << strerror(errno) << "\n";
        }
    }
    

//TILE SIDELENGTH = 20nm
//initialize the cell grid
cell grid[Nx][Ny];

vector<cell> strained_tiles;

//Cycle through every point on the grid (INITIALIZATION LOOP) this loop runs only once
for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
        //Grid points take the form (i,j) and we're assigning numerical values to the ordered pairs stored in the matrix
        grid[i][j].x = i;
        grid[i][j].y = j;
        grid[i][j].strained = false; // Default: unstrained tile
        
        /*
        // Manually set strained regions
        if (i == 20 && j == 40) {
            grid[i][j].strained = true;
            //Set width of potential well
            grid[i][j].sigma_strained = 7;
            //Set depth of the potential well (1 is default)
            grid[i][j].amplitude = 1;

            strained_tiles.push_back(grid[i][j]);

        }
        if (i == 20 && j == 20) {
            grid[i][j].strained = true;
            grid[i][j].sigma_strained = 7;
            grid[i][j].amplitude = 1;
            strained_tiles.push_back(grid[i][j]);
        }

        if (i == 80 && j == 75) {
            grid[i][j].strained = true;
            grid[i][j].sigma_strained = 7;
            grid[i][j].amplitude = 1;
            strained_tiles.push_back(grid[i][j]);
        }
        */

        // Default spawn probability is 0
        double Gauss_height = 0;
        //grab height of gaussian (used to get spawn probability)
        Gauss_height = bivariate_gaussian(grid[i][j].x, grid[i][j].y);
        
        //determine whether or not an exciton spawns based on the statistical distribution (from bivariate gaussian)
        spawn_exciton(Gauss_height, grid[i][j].x, grid[i][j].y);
        
    }   
}

//POTENTIAL ENERGY LANDSCAPE INITILIZATION
//Calculate how nearby strained regions impact excitons on each tile
for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
        //Calculate the force applied on the tile by all strained regions
        for (int q = 0; q < strained_tiles.size(); q++) {
            double xDiff = strained_tiles.at(q).x - grid[i][j].x;
            double yDiff = strained_tiles.at(q).y - grid[i][j].y;
            double radius_strain = sqrt(((xDiff * xDiff) + (yDiff * yDiff)));

            //Calculate force experienced by excitons on the tile by strained_tile(q)

            double force_gaussian_exponent = -( (radius_strain * radius_strain) / (2 * strained_tiles.at(q).sigma_strained * strained_tiles.at(q).sigma_strained)); 

            double Force_Strained_q = K_sA*strained_tiles.at(q).amplitude*(pow(2.71828, force_gaussian_exponent));

            //Calculate components of force
            double F_x = Force_Strained_q*(xDiff/radius_strain);
            double F_y = Force_Strained_q*(yDiff/radius_strain);

            //Add net force and net force components to the struct elements of that tile
            grid[i][j].F_net += Force_Strained_q;
            grid[i][j].F_net_x += F_x;
            grid[i][j].F_net_y += F_y;
        }

        //Calculate the probability of exciton motion from total force on the tile from all strained regions
        double Force_normalizer = abs(grid[i][j].F_net_x) + abs(grid[i][j].F_net_y);

        if (grid[i][j].F_net_x < 0){

                //Calculate probability of exciton moving left on this tile (normalized by denominator)
                grid[i][j].P_left = (1 - grid[i][j].F_net_x)/(4 + Force_normalizer);
                grid[i][j].P_right = 1/(4 + Force_normalizer);
            }

            if (grid[i][j].F_net_x >= 0){

                //Calculate probability of exciton moving left on this tile (normalized by denominator)
                grid[i][j].P_right = (1 + grid[i][j].F_net_x)/(4 + Force_normalizer);
                grid[i][j].P_left = 1/(4 + Force_normalizer);
            }

            if (grid[i][j].F_net_y < 0){

                //Calculate probability of exciton moving down on this tile (normalized by denominator)
                grid[i][j].P_down = (1 - grid[i][j].F_net_y)/(4 + Force_normalizer);
                grid[i][j].P_up = 1/(4 + Force_normalizer);
            }

            if (grid[i][j].F_net_y >= 0){

                //Calculate probability of exciton moving down on this tile (normalized by denominator)
                grid[i][j].P_up = (1 + grid[i][j].F_net_y)/(4 + Force_normalizer);
                grid[i][j].P_down = 1/(4 + Force_normalizer);
            }
    }
}

//define variables for the main loop
int index1 = 1;
double x_expectation_value_sum = 0;
double x_squared_expectation_value_sum = 0;
double y_expectation_value_sum = 0;
double y_squared_expectation_value_sum = 0;
variance variance_total;
double simulation_time = 0;
double converted_variance_average = 0;

//CREATE EXTERNAL FILES THAT
//Define variance file outside the loop (cause we only want 1 variance file)
string file_path_var = variance_path + "/Variance.csv";
std::ofstream varianceFile;
varianceFile.open(file_path_var);
// Write headers for columns
varianceFile << "Variance, Time, \n";


//Define regional density file outside the loop (cause we only want 1 variance file)
string file_path_var2 = regional_density_path1 + "/Exciton_Regional_Density.csv";
std::ofstream RegionalDensityFile1;
RegionalDensityFile1.open(file_path_var2);
// Write headers for columns
RegionalDensityFile1 << "Density, Time, \n";


//get the initial radii of all excitons
for (int i = 0; i < active_Excitons.size(); i++) {
    //Get exciton radii to calculate the variance
    active_Excitons.at(i) = get_radii(active_Excitons.at(i));
}

//MAIN LOOP (590 cycles = 1 nanosecond, 885 cycles = 1.5ns)
int j=0;
while (j<590){

    //Get the total # of living excitons
    int exciton_population_size = active_Excitons.size();

    //Get the variance for this cycle and get the average between x and y variance (since it's ~symmetric)
    variance_total = calculate_variance(variance_total);

    double variance_average = (variance_total.x_variance_struct + variance_total.y_variance_struct)/2;

    if (dimensionless == false){
        //convert variance from units of tiles (:= 400nm^2) to square micrometers
        variance_average *= 0.0004;
    }

    //cout<<"avg variance: "<<variance_average<<"\n";

    //Store the average variance in the variance file every iteration along with the cycle number it was from
    varianceFile << variance_average << ", " << simulation_time << "\n";


    //REGIONAL EXCITON DENSITY CALCULATION (switch on/off -- region1 bool -- see global vars)
        //Customize the region boundaries by editing the exciton_density_region_analysis1 function
    if (region1){
        //Get the density of excitons in units of [excitons/unit tile]
        double Region1_density = exciton_density_region_analysis1();

        //Store the average variance in the variance file every iteration along with the cycle number it was from
        RegionalDensityFile1 << Region1_density << ", " << simulation_time << "\n";
    }

    

    if (j%5 == 0){
        //POSITION DATA STORAGE
        // Create a new csv file for the exciton locations
        string file_path_pos = positons_path + "/Excitons_";
        
        //Add a 1 to the end of the data name for each output (impractical temporary solution to Rstudio bug)
        for (int i=0; i<index1; i++){
            file_path_pos += std::to_string(1);
        }
        file_path_pos += + ".csv"; 
        std::ofstream excitonFile;
        excitonFile.open(file_path_pos);

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
    
        //Set new gaussian variance calculation variables to 0 (allows for new calculation in the main loop)
        x_expectation_value_sum = 0;
        x_squared_expectation_value_sum = 0;
        y_expectation_value_sum = 0;
        y_squared_expectation_value_sum = 0;
        variance_total.x_variance_struct = 0;
        variance_total.y_variance_struct = 0;

    //Cycle through all living excitons
    for (int i = 0; i < exciton_population_size; i++) {
        if (active_Excitons.at(i).active){

            //Determine which way the exciton will move based on probabiilities associated with the tile it's on
                //Instead of grid[i][j] I'm getting the coordinates of the exciton's X and Y positions since those will match that of the tile it's on and putting that into the [i][j] arguments
            int Movement_direction = weightedRandomNum1_4(grid[active_Excitons.at(i).X_pos][active_Excitons.at(i).Y_pos].P_right, 
            grid[active_Excitons.at(i).X_pos][active_Excitons.at(i).Y_pos].P_up, grid[active_Excitons.at(i).X_pos][active_Excitons.at(i).Y_pos].P_left,
            grid[active_Excitons.at(i).X_pos][active_Excitons.at(i).Y_pos].P_down);

            //Move the exciton's position on the grid
            active_Excitons.at(i) = move_exciton(Movement_direction, active_Excitons.at(i));

            //cout<<"Geoge moved to: ("<<active_Excitons.at(i).X_pos<<", "<<active_Excitons.at(i).Y_pos<<")"<<endl;

            //Add 1 to the number of cycles the exciton has been alive
            active_Excitons.at(i).decayTime++;

            //Update exciton radii to calculate the variance
            active_Excitons.at(i) = get_radii(active_Excitons.at(i));

            

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
        //cout<<"Code block run \n";
        //increase while-loop index
        j++;
        //increase time by 1 unit (1 cycle = 0.001695 nanoseconds)
        if (dimensionless == false){
            simulation_time += 0.001695;
        }
        else{
            //Use unit cycles for time in the dimensionless case
            simulation_time += 1;
        }

}

    varianceFile.close();
    RegionalDensityFile1.close();

    if (dimensionless == false){
        cout<<"Simulated "<<simulation_time<<" nanoseconds of diffusion \n";
    }
    else{
        cout<<"Simulated "<<simulation_time<<" simulation cycles \n";
    }
}


//This function generates our gaussian distribution stamp for 1 grid point
double bivariate_gaussian(double x, double y) {
    //Bivariate gaussian parameters
    const double A = 1.0; // amplitude (if A != 1 then statistics will break)
    const double x0 = (Nx-1)/2.0; // center x -- default to halfway along x grid
    const double y0 = (Ny-1)/2.0; // center y -- default to halfway along x grid (this centers the laser impulse at the center of the grid)
    const double sigma_x = 10; // standard deviation along x-axis
    const double sigma_y = 10; // standard deviation along y-axis
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

int weightedRandomNum1_4(double P_right, double P_up, double P_left, double P_down){
    std::random_device rd;
    std::mt19937 rngw(rd());
    std::vector<double> weights;
    weights.clear();
    weights.push_back(P_right);
    weights.push_back(P_up);
    weights.push_back(P_left);
    weights.push_back(P_down);

    std::discrete_distribution<> dist(weights.begin(), weights.end());

    int result = dist(rngw) +1;

    return result;
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


Exciton decay_chance(Exciton exciton){
    //randomly generate a number between 1 and 500 (or whatever the spontaneous decay rate is going to be)
    int randDecayNum;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(1, 500);
    randDecayNum = dist(mt);

    //If we roll a 1, the exciton spontaneously decays
    if (randDecayNum == 1){
        exciton.active = false;
        exciton.decayX = exciton.X_pos;
        exciton.decayY = exciton.Y_pos;
        exciton.fate = "QE before localization";
        
        return exciton;
    }

    else{
        return exciton;
    }
}

Exciton get_radii(Exciton exciton){
    double x = exciton.X_pos;
    double y = exciton.Y_pos;

    //get gaussian center coordinates
    const double x0 = (Nx-1)/2.0;
    const double y0 = (Ny-1)/2.0;

    //Get distance from center along x
    exciton.expectation_x_term1 = (x - x0)*(x - x0);
    exciton.expectation_x_term2 = (x - x0);

    //Get distance from center along y
    exciton.expectation_y_term1 = (y - y0)*(y - y0);
    exciton.expectation_y_term2 = (y - y0);

    // cout<<"x_term1: "<<exciton.expectation_x_term1<<"\n x_term_2: "<<exciton.expectation_x_term2<<"\n";
    // cout<<"y_term1: "<<exciton.expectation_y_term1<<"\n y_term_2: "<<exciton.expectation_y_term2<<"\n";

    return exciton;
}

variance calculate_variance(variance variance_total){
    double sum_x_term1 = 0.0;
    double sum_x_term2 = 0.0;
    double sum_y_term1 = 0.0;
    double sum_y_term2 = 0.0;
    double num_excitons = active_Excitons.size();

    for(int i=0; i<num_excitons; i++){
        Exciton exciton = active_Excitons[i];
        sum_x_term1 += exciton.expectation_x_term1;
        sum_x_term2 += exciton.expectation_x_term2;
        sum_y_term1 += exciton.expectation_y_term1;
        sum_y_term2 += exciton.expectation_y_term2;
    }


    //Calculate variance along x
    double mean_x_term1 = sum_x_term1/num_excitons;
    double mean_x_term2 = sum_x_term2/num_excitons;
    variance_total.x_variance_struct = mean_x_term1 - (mean_x_term2 * mean_x_term2);

    //Calculate variance along y
    double mean_y_term1 = sum_y_term1/num_excitons;
    double mean_y_term2 = sum_y_term2/num_excitons;
    variance_total.y_variance_struct = mean_y_term1 - (mean_y_term2 * mean_y_term2);
    
    return variance_total;
}

double exciton_density_region_analysis1(void){

    //Customize parameters of region
    
    //Set the boundaries of the region (should not exceed bounds of the Nx by Ny region)
    const int x_min = 30;
    const int x_max = 50;
    const int y_min = 30;
    const int y_max = 50;
    
    // Calculate the total area of the region
    double region_width = x_max - x_min;
    double region_height = y_max - y_min;
    double area = region_width * region_height;

    if (area == 0){
        cout<<"ERROR: zero area \n";
        return NAN;
    }

    double num_excitons_in_region = 0;

    //Check all excitons to see which are in the region and then count them
    for (int n=0; n<active_Excitons.size(); n++){
        if (active_Excitons.at(n).X_pos <= x_max){
            if (active_Excitons.at(n).X_pos >= x_min){
                if (active_Excitons.at(n).Y_pos <= y_max){
                    if (active_Excitons.at(n).Y_pos >= y_min){
                        num_excitons_in_region += 1;

                    }
                }
            }
        }
    }

    double density_of_excitons = num_excitons_in_region/area;

    //cout<<"Density = "<<num_excitons_in_region<<"/"<<area<<" = "<<density_of_excitons<<endl;

    return density_of_excitons;

}





/*
double total = grid[25][49].P_left + grid[25][49].P_right + grid[25][49].P_up + grid[25][49].P_down;

cout<<"Tile: ("<<grid[23][41].x<<", "<<grid[23][41].y<<") has F_net of "<<grid[25][49].force<<endl;
cout<<"Tile: ("<<grid[23][41].x<<", "<<grid[23][41].y<<") has F_x of "<<grid[25][49].F_x<<endl;
cout<<"Tile: ("<<grid[23][41].x<<", "<<grid[23][41].y<<") has F_y of "<<grid[25][49].F_y<<endl;

cout<<"Tile: ("<<grid[23][41].x<<", "<<grid[23][41].y<<") has P_left of "<<grid[25][49].P_left<<endl;
cout<<"Tile: ("<<grid[23][41].x<<", "<<grid[23][41].y<<") has P_right of "<<grid[25][49].P_right<<endl;
cout<<"Tile: ("<<grid[23][41].x<<", "<<grid[23][41].y<<") has P_up of "<<grid[25][49].P_up<<endl;
cout<<"Tile: ("<<grid[23][41].x<<", "<<grid[23][41].y<<") has P_down of "<<grid[25][49].P_down<<endl;


//So clearly there is an issue regarding the way I've gotten the force components
cout<<"Sum of probabilities: "<<total<<endl;

*/


/*
for (int i=0; i<20; i++){
//RANDOM EXCITON MOTION

cout<<"("<<active_Excitons.at(0).X_pos<<", "<<active_Excitons.at(0).Y_pos<<") \n";
int random_move_num = randomNum1_4();
active_Excitons.at(0) = move_exciton(random_move_num, active_Excitons.at(0));
}
*/


//SINGLE TEST EXCITON
//Instantiate test exciton(s) [I put it in this so that the exciton is only instantiated once, not a very clean way to code, David]
//Exciton test_exciton1(20,30,true,"George");

//Add the exciton to the list of active (not decayed) excitons
//active_Excitons.push_back(test_exciton1);
