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
    bool strain_dectection;

    Exciton(int X_pos, int Y_pos, bool active, std::string name, bool strain_dectection){
        this->X_pos = X_pos;
        this->Y_pos = Y_pos;
        this->active = active;
        this->name = name;
        decayTime=0;
        this->strain_dectection = strain_dectection;
    }
};

//GLOBAL
using std::cout;
using std::endl;
using std::string;
using std::vector;
int spawn_index = 0;

//Strain attraction parameter
const double K_sA = 500;

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
    bool strained;
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
Exciton get_radii(Exciton exciton);
variance calculate_variance(variance variance_total);
int weightedRandomNum1_4(double P_right, double P_up, double P_left, double P_down);


int main(){
    //CREATE FOLDERS FOR THE SIMULATION OUTPUT
    // Get current date and time
    std::time_t t = std::time(nullptr);
    char date_time[100];
    std::strftime(date_time, sizeof(date_time), "%Y-%m-%d %H:%M:%S", std::localtime(&t));

    // Create directory path with unique name
    string positons_path = "/Users/david/Desktop/C++/Exciton Research/Heat Map R/" + string(date_time) + string(" Positons");
    string variance_path = "/Users/david/Desktop/C++/Exciton Research/Heat Map R/" + string(date_time) + string(" Variance");

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

//TILE SIDELENGTH = 20nm
//initialize the cell grid
cell grid[Nx][Ny];

//Cycle through every point on the grid (INITIALIZATION LOOP) this loop runs only once
for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
        //Grid points take the form (i,j) and we're assigning numerical values to the ordered pairs stored in the matrix
        grid[i][j].x = i;
        grid[i][j].y = j;
        grid[i][j].strained = false; // Default: unstrained tile
        
        // Manually set strained regions
        if (i == 20 && j == 20) {
            grid[i][j].strained = true;
        }
        if (i == 80 && j == 20) {
            grid[i][j].strained = true;
        }
        if (i == 20 && j == 80) {
            grid[i][j].strained = true;
        }
        if (i == 80 && j == 80) {
            grid[i][j].strained = true;
        }
        
        
    

        // Default spawn probability is 0
        double Gauss_height = 0;
        //grab height of gaussian (used to get spawn probability)
        Gauss_height = bivariate_gaussian(grid[i][j].x, grid[i][j].y);
        
        //determine whether or not an exciton spawns based on the statistical distribution (from bivariate gaussian)
        spawn_exciton(Gauss_height, grid[i][j].x, grid[i][j].y);
        
        //Instantiate test excitons
        //Exciton test_exciton1(40,40,true,"george",true);

        //Add the exciton to the list of active (not decayed) excitons
        //active_Excitons.push_back(test_exciton1);
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

double F_right, F_left, F_up, F_down, F_net;
double P_right, P_left, P_up, P_down;

//Define cell struct for the strained tile
cell detected_strainCell;
vector<cell> strained_tiles;

//Define variance file outside the loop (cause we only want 1 variance file)
string file_path_var = variance_path + "/Variance.csv";
std::ofstream varianceFile;
varianceFile.open(file_path_var);
// Write headers for columns
varianceFile << "Variance, Time, \n";

//get the initial radii of all excitons for the t=0 variance calculation
for (int i = 0; i < active_Excitons.size(); i++) {
    //Get exciton radii to calculate the variance
    active_Excitons.at(i) = get_radii(active_Excitons.at(i));
}

    //MAIN LOOP (590 cycles = 1 nanosecond, 885 cycles = 1.5ns)
    int j=0;
    while (j<500){
        //Get the total # of living excitons
        int exciton_population_size = active_Excitons.size();

        //Get the variance for this cycle and get the average between x and y variance (since it's ~symmetric)
        variance_total = calculate_variance(variance_total);

        double variance_average = (variance_total.x_variance_struct + variance_total.y_variance_struct)/2;

        //convert variance from units of tiles (:= 400nm^2) to square micrometers
        converted_variance_average = variance_average*0.0004;

        //cout<<"avg variance: "<<variance_average<<"\n";

        //Store the average variance in the variance file every iteration along with the cycle number it was from
        varianceFile << converted_variance_average << ", " << simulation_time << "\n";
       
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
//TURNING OFF RANDOM DECAY AND ANNIHILATION 
/*
            //Chance to decay randomly
            active_Excitons.at(i) = decay_chance(active_Excitons.at(i));
            
            //Recombination check (destroy two active excitons if they share the same tile)
            for (int k = 0; k < exciton_population_size; k++) {
                if (active_Excitons.at(i).active && active_Excitons.at(k).active){
                    if(active_Excitons.at(i).name != active_Excitons.at(k).name){
                        if (active_Excitons.at(i).X_pos == active_Excitons.at(k).X_pos && active_Excitons.at(i).Y_pos == active_Excitons.at(k).Y_pos){
                            //cout<<"Destroyed "<< active_Excitons.at(i).name<<" and "<<active_Excitons.at(k).name<<"\n";
                            //cout<<active_Excitons.at(i).name<<" coordinates are ("<<active_Excitons.at(i).X_pos<<","<<active_Excitons.at(i).Y_pos<<") and "<<active_Excitons.at(k).name<<" coordinates are ("<<active_Excitons.at(k).X_pos<<","<<active_Excitons.at(k).Y_pos<<") \n";
                            active_Excitons.at(i).active = false;
                            active_Excitons.at(k).active = false;

                            //Store time of death stats
                            active_Excitons.at(i).decayX = active_Excitons.at(i).X_pos;
                            active_Excitons.at(i).decayY = active_Excitons.at(i).Y_pos;
                            active_Excitons.at(i).fate = "Recombination";
                            active_Excitons.at(k).decayX = active_Excitons.at(k).X_pos;
                            active_Excitons.at(k).decayY = active_Excitons.at(k).Y_pos;
                            active_Excitons.at(k).fate = "Recombination";
                        }
                    }
                }
            }
            */
        //REMOVE TRAPPED EXCITONS FROM THE SYSTEM
        // if (active_Excitons.at(i).X_pos == 20 && active_Excitons.at(i).Y_pos == 80){
        //     active_Excitons.at(i).active = false;
        //     active_Excitons.at(i).fate = "Trapped in strain";

        // }
        // if (active_Excitons.at(i).X_pos == 80 && active_Excitons.at(i).Y_pos == 20){
        //     active_Excitons.at(i).active = false;
        //     active_Excitons.at(i).fate = "Trapped in strain";
        // }
        // if (active_Excitons.at(i).X_pos == 80 && active_Excitons.at(i).Y_pos == 80){
        //     active_Excitons.at(i).active = false;
        //     active_Excitons.at(i).fate = "Trapped in strain";
        // }
        // if (active_Excitons.at(i).X_pos == 20 && active_Excitons.at(i).Y_pos == 20){
        //     active_Excitons.at(i).active = false;
        //     active_Excitons.at(i).fate = "Trapped in strain";
        // }
        
        //STRAINED EVOLUTION
        if (active_Excitons.at(i).active){
            //Search for nearby strained regions (default 1 tile=20nm)
            //Get the exciton's cell location
            int centerSearch_X = active_Excitons.at(i).X_pos;
            int centerSearch_Y = active_Excitons.at(i).Y_pos;
            //Search the nearest 10 tiles in all directions
            for (int x = centerSearch_X - 30; x <= centerSearch_X + 30; x++) {
                for (int y = centerSearch_Y - 30; y <= centerSearch_Y + 30; y++) {
                    // Check if the current tile is within the bounds of the grid
                    if (x >= 0 && x < Nx && y >= 0 && y < Ny) {
                        // Check if the current tile is strained
                        if (grid[x][y].strained == true) {
                            //cout << "Detected nearby strain at (" << x << "," << y << ")\n";

                            //identify the strained tile
                            detected_strainCell.x = grid[x][y].x;
                            detected_strainCell.y = grid[x][y].y;
                            detected_strainCell.strained = true;
                            
                            // Add the strained tile to a list or perform some other action
                            strained_tiles.push_back(detected_strainCell);
                        }
                    }
                }
            }
            //cout<<"block 1 done \n";
            if (strained_tiles.size() == 0){
                active_Excitons.at(i).strain_dectection = false;
            }
            else{
                //Reset forces to 0 for a new calculation
                F_down = 0;
                F_up = 0;
                F_right = 0;
                F_left = 0;
                for (int q = 0; q < strained_tiles.size(); q++) {
                    double xDiff = strained_tiles.at(q).x - active_Excitons.at(i).X_pos;
                    double yDiff = strained_tiles.at(q).y - active_Excitons.at(i).Y_pos;
                    double radius_strain = sqrt(((xDiff * xDiff) + (yDiff * yDiff)));
                    
                    //calculate "force" along the up,down,left,right axes
                    if (xDiff < 0) {
                        //multiply by -1 to make F positive
                        F_left += (K_sA * xDiff * -1)/(radius_strain * radius_strain * radius_strain);            
                    } 
                    else if (xDiff > 0){
                        F_right += (K_sA * xDiff)/(radius_strain * radius_strain * radius_strain);
                    }

                    if (yDiff < 0) {
                        //multiply by -1 to make F positive
                        F_down += (K_sA * yDiff * -1)/(radius_strain * radius_strain * radius_strain);
                    } 
                    else if (yDiff > 0){
                        F_up += (K_sA * yDiff)/(radius_strain * radius_strain * radius_strain);
                    }
                }   
            
            F_net = F_right + F_down + F_left + F_up;
            //calculate probability of moving in a direction
            P_right = (1 + F_right)/(4 + F_net);
            P_left = (1 + F_left)/(4 + F_net);
            P_up = (1 + F_up)/(4 + F_net);
            P_down = (1 + F_down)/(4 + F_net);

            //cout<<"P_r= "<<P_right<<", P_u= "<<P_up<<"P_l= "<<P_left<<", P_d= "<<P_down<<"\n";

            //Randomly generate with integer from 1 to 4 with those weighted probabilities
            int random_strained_direction;
            random_strained_direction = weightedRandomNum1_4(P_right, P_up, P_left, P_down);

            //Move the exciton
            active_Excitons.at(i) = move_exciton(random_strained_direction, active_Excitons.at(i));

            //Add 1 to the number of cycles the exciton has been alive
            active_Excitons.at(i).decayTime++;

            }
        }
        //Reset detected tiles for the next exciton
        strained_tiles.clear();
        
        //RANDOM EVOLUTION (for unstrained or no nearby strain detected)
        if (active_Excitons.at(i).strain_dectection == false){

            //Check to see if the exciton has decayed, if so, don't run the time evolution code block 
            if (active_Excitons.at(i).active){
                //Generate a random number between 1 and 4
                int randomDirection;
                randomDirection = randomNum1_4();

                //Move the exciton in a totally random direction
                active_Excitons.at(i) = move_exciton(randomDirection, active_Excitons.at(i));

                //Add 1 to the number of cycles the exciton has been alive
                active_Excitons.at(i).decayTime++;

                //cout<<"Moved randomly"<<endl;
            }
            
        }
        //Update exciton radii to calculate the variance
        active_Excitons.at(i) = get_radii(active_Excitons.at(i));

        //reset strain detection (allows us to go do a new scan every cycle)
        active_Excitons.at(i).strain_dectection = true;
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
        simulation_time += 0.001695;
    }
    varianceFile.close();
    cout<<"Simulated "<<simulation_time<<" nanoseconds of diffusion \n";
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
        Exciton exciton(x,y,true,exciton_name,true);

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

