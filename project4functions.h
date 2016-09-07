
#include <math.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <deque>

using namespace std;

struct Coordinate{

	int x;
	int y;
	char type; //'w' = wild, 'b' = border, 's' =  safe 

};

struct OptCoord{
  int x;
  int y; 
  bool visited; 
};

struct Prims{
  bool visited; //whether or not it has been visited
  double distance; // distance 
  int parent; //index of parent 
};





double mst_distance(Coordinate &a, Coordinate &b)
{
	if (a.type == 'w' && b.type == 's') //one is safe and the other is wild 
	{
		return std::numeric_limits<double>::infinity();
	}
	else if (b.type == 'w' && a.type == 's')
	{
		return std::numeric_limits<double>::infinity(); 
	}

	return sqrt(1.0*(b.x-a.x)*(b.x-a.x)+1.0*(b.y-a.y)*(b.y-a.y));
}

double euclidean_distance(OptCoord &a, OptCoord &b)
{
  return sqrt(1.0*(b.x-a.x)*(b.x-a.x)+1.0*(b.y-a.y)*(b.y-a.y));
}


double prims_algorithm(ostream &ss, int path_size, vector<Coordinate> &coordinates, vector<Prims> &prims_guys)
{
  
  prims_guys[0].distance = 0; 

  double min_distance = std::numeric_limits<double>::infinity(); 
  double total_distance = 0; 
  bool all_visited = false; 
  int visited_counter = 0; 
  int min_index = 0;  

  while (all_visited != true)
  {   
 	
    for(int i = 0; i < path_size; i++) //find current shortest distance
    {
      if(prims_guys[i].visited == false)
      {
  			if (prims_guys[i].distance < min_distance) 
  			{
  				min_distance = prims_guys[i].distance;
  				min_index = i; 	
  			} 
      }   
    }//end of for loop 

    prims_guys[min_index].visited = true; 
    visited_counter++; 

    if(visited_counter == path_size)
    {
    	all_visited = true; 
    }

    if(min_index < prims_guys[min_index].parent && min_index != 0) //&& all_visited == false)
    {
    	  ss << min_index << " " << prims_guys[min_index].parent << '\n';
    }
    else if (min_index != 0)// (all_visited == false)
    {
        ss << prims_guys[min_index].parent << " " << min_index << '\n';
    }
    
    for(int j = 0; j < path_size; j++) //loop thru 
    {
    	if(prims_guys[j].visited == false) //&& can_connect(coordinates[min_index], coordinates[k])) ///
    	{
    		if(mst_distance(coordinates[j], coordinates[min_index]) < prims_guys[j].distance)
    		{
    			prims_guys[j].distance = mst_distance(coordinates[j], coordinates[min_index]);
    			prims_guys[j].parent = min_index; 
    		}
    	}
    }   
    total_distance = total_distance + min_distance; 
   
   min_distance = std::numeric_limits<double>::infinity(); 

 }//end of while 

 return total_distance; 



}


double twoOpt(double total_distance, vector<OptCoord> &coords, vector<int> &neighbors)
{

  vector<int>::iterator begin_reverse; //= neighbors.begin();
  vector<int>::iterator end_reverse;// = neighbors.begin(); 

  //cout << "neighbors size " << (int)neighbors.size() << endl; 

  for(int i = 0; i < (int)neighbors.size(); i++)
  {
    for(int j = i + 2; j < (int)neighbors.size() - 1; j++) //CHECK THIS!!!!!!!
    {
      
      if((euclidean_distance(coords[neighbors[i]], coords[neighbors[i+1]]) 
        + euclidean_distance(coords[neighbors[j]], coords[neighbors[j+1]])) >
        (euclidean_distance(coords[neighbors[i]], coords[neighbors[j]]) 
         + euclidean_distance(coords[neighbors[i+1]], coords[neighbors[j+1]])))
      { 
        
          //subtract old sum 
          total_distance = total_distance - (euclidean_distance(coords[neighbors[i]], coords[neighbors[i+1]]) 
          + euclidean_distance(coords[neighbors[j]], coords[neighbors[j+1]]));
          
          //add new sum 
          total_distance = total_distance + (euclidean_distance(coords[neighbors[i]], coords[neighbors[j]]) 
           + euclidean_distance(coords[neighbors[i+1]], coords[neighbors[j+1]])); 

          begin_reverse = neighbors.begin() + i + 1;
          end_reverse = neighbors.begin() + j + 1; 
          reverse(begin_reverse, end_reverse);
        //add new sum 
      }
    }

  }

  return total_distance; 

}


double nearest_neighbor(int path_size, vector<OptCoord> &coordinates, vector<int> &neighbors)
{

  double min_distance = std::numeric_limits<double>::infinity();
  int min_index = 0; 
  coordinates[0].visited = true; //mark first index visited
  int visited = 1; //initial point is visited 
  bool all_visited = false; 
  int previous_coordinate = 0; 
  double total_distance = 0; 

  int counter = 0; 
  
  int first_index = 0;

  for(int i = 1; i < (int)coordinates.size(); i++) //find the distance between 0 and closest point 
  {
    if(euclidean_distance(coordinates[0], coordinates[i]) < min_distance)
    {
      min_distance = euclidean_distance(coordinates[0], coordinates[i]);
      first_index = i; //assign first index to first closest point 
    }
  }

  total_distance = min_distance; //add to total distance this initial distance 
  coordinates[first_index].visited = true; //mark that index visited
  visited++; //2 points are visitied 

  //push both back into the vector 
  neighbors.push_back(0);
  neighbors.push_back(first_index);

  previous_coordinate = first_index; 
  min_distance = std::numeric_limits<double>::infinity(); //reset min distance 

  while(!all_visited)
  {   
      
      for(int i = 0; i < (int)coordinates.size(); i++)
      {
     
        if(coordinates[i].visited == false)
        {
         
          if(euclidean_distance(coordinates[i], coordinates[previous_coordinate]) < min_distance)
          {

            min_distance = euclidean_distance(coordinates[i], coordinates[previous_coordinate]);

            min_index = i; //find the next neighbor point that is closest
            
          }
        }
     }  


  neighbors.push_back(min_index); //push that point back 
  coordinates[min_index].visited = true; //mark it visited 
  visited++; 
  previous_coordinate = min_index; //that point becomes the next ~previous coordinate~
  total_distance = total_distance + min_distance; //update TOTAL DISTANCE WITH NEW ONE 
 
  counter++; 
  min_distance = std::numeric_limits<double>::infinity(); //reset min distance 

    if(visited == path_size)
    {
      all_visited = true; 
    }
  }//end of while loop 

  //add in the distance between the last coordinate and the beginning 
  total_distance = total_distance + euclidean_distance(coordinates[0], coordinates[previous_coordinate]);
  neighbors.push_back(0); //add a zero on to the end of the path 

  
  total_distance = twoOpt(total_distance, coordinates, neighbors);
  
  return total_distance; 
}



//OPTTSP

double mst(vector<vector<double>> &dist_matrix, vector<double> &mst_distances, 
vector<OptCoord> &coordinates, deque<int> &unvisited)
{
 // cout << "ENTERING MST " << endl; 
  for(int i = 0; i < (int)coordinates.size(); i++) //set all of them back to false
   {
     coordinates[i].visited = false; 
     mst_distances[i] = std::numeric_limits<double>::infinity();
   }
   /*
  cout << "UNVISITED SIZE " << unvisited.size() << endl;
  cout << "printing UNVISTED IN MST " << endl;
  for(int i = 0; i < (int)unvisited.size(); i++)
  {
    cout << unvisited[i] << " ";
  }
  cout << endl; 
  cout << " ^PRINTED^ " << endl; 
 */
  
  mst_distances[unvisited[0]] = 0; 
  double min_distance = std::numeric_limits<double>::infinity(); 
  double total_distance = 0; 
  bool all_visited = false; 
  int visited_counter = 0; 
  int min_index = 0;  

  while (all_visited != true)
  {    
    //cout << (int)unvisited.size() << endl;
  
    for(int i = 0; i < (int)unvisited.size(); i++) //find current shortest distance
    {
      if(coordinates[unvisited[i]].visited == false)
      {
       // cout << mst_distances[unvisited[i]] << endl; 
        if (mst_distances[unvisited[i]] < min_distance) 
        {
          min_distance = mst_distances[unvisited[i]];
          min_index = unvisited[i];  
        } 
      }   
    }//end of for loop 
    
    //cout << "here " << endl;
    //cout << min_index << endl; 
    coordinates[min_index].visited = true; 
    visited_counter++; 

    if(visited_counter == (int)unvisited.size())
    {
      all_visited = true; 
    }
    
    for(int j = 0; j < (int)unvisited.size(); j++) //loop thru 
    {
      if(coordinates[unvisited[j]].visited == false)
      {
        if(dist_matrix[unvisited[j]][min_index] < mst_distances[unvisited[j]])
        {
          mst_distances[unvisited[j]] = dist_matrix[unvisited[j]][min_index];
        }
      }
    } 
  //cout << "min distance " << min_distance << endl; 
  //cout << "total dist " << total_distance << endl; 
  total_distance = total_distance + min_distance; 
  //cout << "total distance added " << total_distance << endl;
  min_distance = std::numeric_limits<double>::infinity(); 

 
 }//end of while 

// cout << "EXITING MST " <<endl; 

 /*
 for(int i = 0; i < (int)coordinates.size(); i++) //set all of them back to false
 {
   coordinates[i].visited = false; 
   mst_distances[i] = std::numeric_limits<double>::infinity();
 }
  */
 return total_distance; 
}

bool promising(vector<double> mst_distances, deque<int> &unvisited, vector<int> &path, 
vector<vector<double>> &dist_matrix, vector<OptCoord> &coordinates, double &best_path_length, 
double current_path_length)
{
  double potential_path_length = 0; 
  

 // cout << "****unvisited: ";
  /*
  for(int i = 0; i < (int)unvisited.size();i++)
  {

    cout << unvisited[i] << " ";
  }
  cout << endl;
   cout << "****path: ";
  for(int i = 0; i < (int)path.size();i++)
  {

    cout << path[i] << " ";
  }
  cout << endl;
*/
  double mst_weight = mst(dist_matrix, mst_distances, coordinates, unvisited);


  //cout << "mst weight " <<  mst_weight << endl; 
  double arm_one = 0; 
  double arm_two = 0;
  double min_distance_begin = std::numeric_limits<double>::infinity();
  double min_distance_end = std::numeric_limits<double>::infinity();

  
  int front_path = path.front();
  int back_path = path.back(); 
  //cout << "PATH FRONT " << front_path << endl; 
 // cout << "PATH BACK " << back_path << endl; 
  
  for(int i = 0; i < (int)unvisited.size(); i++)
  {
     

    if(dist_matrix[front_path][unvisited[i]] < min_distance_begin)
    {
      //cout << "****DISTANCE FRONT point " << unvisited[i] << " and " << front_path << ": " << dist_matrix[front_path][unvisited[i]] << endl;
      min_distance_begin = dist_matrix[front_path][unvisited[i]];
      //cout << min_distance_begin << endl;
    }
    if(dist_matrix[back_path][unvisited[i]] < min_distance_end)
    {
      //cout << "****DISTANCE BACK point " << unvisited[i] << " and " << back_path << ": " << dist_matrix[back_path][unvisited[i]] << endl;
    
      
      min_distance_end = dist_matrix[back_path][unvisited[i]];
      //cout << min_distance_end << endl; 
    }
  }

  arm_one = min_distance_begin;
  arm_two = min_distance_end;
  //cout << "arm one " << arm_one << endl; 
  //cout << "arm two " << arm_two << endl;


  //cout << "current path length " << current_path_length << endl; 
  potential_path_length = current_path_length + arm_one + arm_two + mst_weight; 
 // cout << "potential_path_length " << potential_path_length << endl;

  if(potential_path_length < best_path_length)
  {
    return true;
  }

  return false; 
}



void genPerms(vector<double> &mst_distances, deque<int> &unvisited, vector<int> &path, vector<int> &best_path, vector<vector<double>> &dist_matrix, 
vector<OptCoord> &coordinates, double &best_path_length, double current_path_length)
{
  if(unvisited.empty())
  {
    current_path_length = current_path_length + dist_matrix[path.front()][path.back()];
    if(current_path_length < best_path_length)
    {
      best_path_length = current_path_length; 
      best_path = path;
    }
    current_path_length = current_path_length - dist_matrix[path.front()][path.back()];
    return; 
  }
  if(!promising(mst_distances, unvisited, path, dist_matrix, coordinates, best_path_length, current_path_length))
  {
    //cout << " not promising " << endl;
    return; 
  }
  for (unsigned k = 0; k != unvisited.size(); k++)
  {
    double new_edge = dist_matrix[unvisited.front()][path.back()];
    //current_path_length = current_path_length + dist_matrix[unvisited.front()][path.back()];
    path.push_back(unvisited.front());
    unvisited.pop_front();
    genPerms(mst_distances, unvisited, path, best_path, dist_matrix, coordinates, best_path_length, current_path_length+new_edge); //??
    unvisited.push_back(path.back());
    path.pop_back();
    //current_path_length = current_path_length - dist_matrix[path.back()][unvisited.back()];
    
    
  }
}













































