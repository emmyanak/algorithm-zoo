
#include <getopt.h>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <limits>
#include <sstream>
#include <iomanip>
#include "project4functions.h"


using namespace std;

int main(int argc, char **argv){

stringstream os; 
os << std::setprecision(2); //show 2 dec places
os << std::fixed; 

cout << std::setprecision(2); //show 2 dec places
cout << std::fixed; 


  static struct option longopts[] = 
    {
       {"mode", required_argument, NULL, 'm'},
       {"help", no_argument, NULL, 'h'},
       {NULL, 0, NULL, 0}
    };

    
    int c;
    int indx = 0;
    char mode = '\0';
 
    while ((c = getopt_long(argc, argv, "m:h", longopts, &indx)) != -1)
    {
        switch(c) { 
            case 'm':
                if (optarg == NULL)
                {
                    cerr << "must provide mode " << '\n';
                    exit(1);
                }
                else
                {
                	if (strcmp(optarg, "MST") == 0)
                	{
               		    mode = 'a';
               	    }
                    else if (strcmp(optarg, "FASTTSP") == 0)
                    {
                   	    mode = 'b';
                    }
                    else if(strcmp(optarg, "OPTTSP") == 0)
                    {
                   	    mode = 'c';
                    }
                    else 
                    {
                    	cerr << "incorrect output mode " << '\n';
                    	exit(1);
                    }    
                }
               break; 
            case 'h': 
            cout << "help message " << '\n'; 
            exit(0);    
        }
    }//end of getopt long while loop 

    

    //read in input (for all)
    vector<Coordinate> zoo_map;
    Coordinate input_coordinate; 	
    vector<OptCoord> opt_input; 
    
    //for MST implementation
    vector<Prims> prims_vector;

    //for FASTOPT & OPTOPT implementation 
    vector<int> neighborz; 
    OptCoord opt_coord; 

    //for OPTOPT
    vector<vector<double>> distance_matrix; 
    vector<double> mst_distances; 
    vector<int> opt_path;
    vector<int> best_path; 
    deque<int> unvisited_nodes;  
     
    Prims mst_point; 

    int num_entries; 
    int first_coord = 0;
    int second_coord = 0;  
    bool border_check = false; 
    bool wild_check = false; 
    bool safe_check = false; 
    
    
    cin >> num_entries;
    //set vector 
     

  if(mode == 'c') //create a distance matrix of size num_entries*num_entries
  {
    //cout << "mode = c" << endl; 
    distance_matrix.resize(num_entries, vector<double>(num_entries));
    //cout << "reset size " << endl; 
  }

    while (cin >> first_coord)
    {

    
        //cout << "first coord " << first_coord << endl; 
    		cin >> second_coord; 

    		if (mode == 'a')//prims algorithm	
    		{
    		    input_coordinate.x = first_coord;
    		    input_coordinate.y = second_coord;
          if(first_coord < 0 && second_coord < 0)
          {
            input_coordinate.type = 'w';//wild point 
            wild_check = true; 
          }
	    		else if((first_coord == 0 && second_coord <= 0) || (first_coord <= 0 && second_coord == 0))
	    		{
	    			input_coordinate.type = 'b'; //border point
	    			border_check = true; 
	    		}
	    		
	    		else
	    		{
	    			input_coordinate.type = 's'; //safe point
            safe_check = true; 
	    		}

	    	  zoo_map.push_back(input_coordinate);

              mst_point.visited = false; 
              mst_point.distance = std::numeric_limits<double>::infinity(); //set distance to infinity 
              mst_point.parent = 0; //set parent to zero 
              prims_vector.push_back(mst_point);
    	    }
    		//if mode == b or c
        
    		opt_coord.x = first_coord;
    		opt_coord.y = second_coord; 
    		opt_coord.visited = false; 


        opt_input.push_back(opt_coord);


        if(mode == 'c')
        {
          mst_distances.push_back(std::numeric_limits<double>::infinity());
          
          /*
          for(int i = 0; i < (int)opt_input.size(); i++) //fill the distance matrix
          {
            distance_matrix[i][(int)opt_input.size()-1] = euclidean_distance(opt_coord, opt_input[i]);
            distance_matrix[(int)opt_input.size()-1][i] = euclidean_distance(opt_input[i], opt_coord);
          }
          */
         
        }
    }
 

   //cout << num_entries << endl; 
   //cout << "mode: " << mode << endl; 


   if((mode == 'a') && (wild_check == true && safe_check == true && border_check == false))
   {
   	cerr << "Cannot construct MST " << endl; 
   	exit(1);
   }
//printint input vector 
   /*
   for(int i = 0; i < (int)opt_input.size(); i++)
   {
	   	cout << "map element " << i << endl; 
	   	cout << opt_input[i].x << ", " << opt_input[i].y << endl; 
	   	cout << opt_input[i].visited << endl; 
   }
  */



   double mst_weight = 0; 
   //cout << "mode " << mode << endl; 

   if(mode == 'a')
   { 
     mst_weight = prims_algorithm(os, num_entries, zoo_map, prims_vector);
     cout << mst_weight << '\n';
     cout << os.str();
   }

   if(mode == 'b')
   {
     
    // cout << "num entries " << num_entries << endl; 
     //cout << "opt input size " << opt_input.size() << endl; 
     //cout << "neighbors size " << neighborz.size() << endl; 
     mst_weight = nearest_neighbor(num_entries, opt_input, neighborz);
     cout << mst_weight << '\n';

     cout << neighborz[0];
     for(int i = 1; i < (int)neighborz.size() - 1; i++)
     {
      cout << " " << neighborz[i]; 
     }

     cout << endl; 
   }

   if(mode == 'c')
   {

    for(int i = 0; i < num_entries; i++) //fill the distance matrix
    {
      for(int j = 0; j < num_entries; j++)
        {      
          distance_matrix[i][j] = euclidean_distance(opt_input[i], opt_input[j]);
          distance_matrix[j][i] = euclidean_distance(opt_input[i], opt_input[j]);
        }
    }

    //test MST 
    /*
    for(int i = 0; i < num_entries; i++)
    {
      unvisited_nodes.push_back(i);
    }

    double test_value = mst(distance_matrix, mst_distances, opt_input, unvisited_nodes);
    
    cout << test_value << endl; 
    */
    

    /*
    for (int i = 0; i < num_entries; i++)
    {
      for(int j = 0; j < num_entries; j++)
      {

        cout << distance_matrix[i][j] << " "; 
      }
      cout << endl; 
    }
    */
    //unvisited_nodes.clear();
   

    
    double best_path_length = nearest_neighbor(num_entries, opt_input, neighborz);
     //cout << best_path_length << endl; 
     //cout << neighborz[0];
     /*
     for(int i = 1; i < (int)neighborz.size() - 1; i++)
     {
      cout << " " << neighborz[i]; 
     }
     */
    double current_path_length = 0; //std::numeric_limits<double>::infinity(); 
    opt_path.push_back(0);
    best_path = neighborz;
    for(int i = 1; i < num_entries; i++)
    {
      unvisited_nodes.push_back(neighborz[i]);
    }
    genPerms(mst_distances, unvisited_nodes, opt_path, best_path, distance_matrix, 
    opt_input, best_path_length, current_path_length);

    cout << best_path_length << endl; 
    cout << best_path[0]; 
    int i_limit = (int)best_path.size();
    if(best_path == neighborz)
    {
      i_limit = (int)best_path.size() - 1; 
    }
    for(int i = 1; i < i_limit; i++)
    {
      cout << " " << best_path[i];
    }
    cout << endl; 
    
   }
   
   





   














}//end of main