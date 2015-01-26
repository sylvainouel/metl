#include <string>
#include <utility>
#include <fstream>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "tsp_prob.hh"

using namespace std;

//this function reads a file
//and compute the distance matrix
//between every city of the problem
tsp_prob::tsp_prob()
  : d(0), prob_size(0), maxrow(0),
    cities(),
    weight_type(TYPE_NONE), 
    candidate(*this)
{}


void tsp_prob::load(const string& tspfile)
{
  ifstream in_file(tspfile.c_str());
  string in_line;
  bool done=false;
  t_weight_format weight_format=FORMAT_NONE;
  unsigned i,j;

  if (in_file.fail()) {
    cerr << "Failed to open file" << endl;
    abort();
    //throw("BAD_FILE");
  }

  while(!done && !in_file.eof()) {
    getline(in_file, in_line, '\n');

    if (in_line.find("DIMENSION")!=string::npos) {
      prob_size = atoi(in_line.substr(in_line.find(":")+1).c_str());
      // << prob_size << endl;
      continue;
    }
    if (in_line.find("TYPE")==0) {
      if (in_line.substr(in_line.find(":")+1).find("TSP")==string::npos)
	{
	  cerr << "Problem type not supported" << endl;
	  in_file.close();
	  abort();
	  //	  throw("BAD_TYPE");
	}
      continue;
    }
    if (in_line.find("NAME")==0) {
      cout << in_line << endl;
      continue;
    }
    if (in_line.find("COMMENT")==0) {
      cout << in_line << endl;
      continue;
    }
    if (in_line.find("EDGE_WEIGHT_TYPE")!=string::npos) {
      if (in_line.substr(18).find("EUC_2D")!=string::npos) {
	weight_type=EUC_2D;
	//cout << "weight type is: EUC_2D" << endl;
	continue;
      }
      if (in_line.substr(18).find("CEIL_2D")!=string::npos) {
	weight_type=CEIL_2D;
	//cout << "weight type is: CEIL_2D" << endl;
	continue;
      }

      if (in_line.substr(18).find("EXPLICIT")!=string::npos) {
	weight_type=EXPLICIT;
	//cout << "weight type is: EXPLICIT" << endl;
	continue;
      }
      if (in_line.substr(18).find("ATT")!=string::npos) {
	weight_type=ATT;
	continue;
      }

    }
    if (in_line.find("EDGE_WEIGHT_FORMAT")!=string::npos) {
      if (in_line.substr(18).find("FULL_MATRIX")!=string::npos) {
	weight_format=FULL_MATRIX;
	//cout << "weight format is: FULL_MATRIX" << endl;
	continue;
      }
      if (in_line.substr(18).find("LOWER_DIAG_ROW")!=string::npos) {
	weight_format=LOWER_DIAG_ROW;
	//cout << "weight format is: LOWER_DIAG_ROW" << endl;
	continue;
      }
      if (in_line.substr(18).find("UPPER_DIAG_ROW")!=string::npos) {
	weight_format=UPPER_DIAG_ROW;
	//cout << "weight format is: UPPER_DIAG_ROW" << endl;
	continue;
      }
      if (in_line.substr(18).find("UPPER_ROW")!=string::npos) {
	weight_format=UPPER_ROW;
	//cout << "weight format is: UPPER_ROW" << endl;
	continue;
      }
    }
    if (in_line.find("NODE_COORD_SECTION")!=string::npos) {
      done=true;
      break;
    }
    if (in_line.find("EDGE_WEIGHT_SECTION")!=string::npos) {
      done=true;
      break;
    }
  }

  if (!done && in_file.eof()) {
    cerr << "bad file" << endl;
    in_file.close();
    abort();
    //    throw("BAD_FILE");
  }

  maxrow = prob_size;
  const unsigned mms = MAX_MATRIX_SIZE==-1 ? 0 : MAX_MATRIX_SIZE; 
  if (MAX_MATRIX_SIZE!=-1 && weight_type!=EXPLICIT)
    {
      maxrow = std::min((unsigned int) prob_size, mms / (unsigned int) (prob_size*sizeof(int)));
    }

  d = new int*[prob_size];
  for (i=0; i<prob_size; i++) {
    if (i>=maxrow) {
      d[i]=0;
    } else {
      d[i] = new int[prob_size];
    }
  }

  switch (weight_type) {
  case TYPE_NONE:
    {
      cerr << "Invalid weight type" << endl;
      in_file.close();
      abort();
      //      throw("INVALID_WEIGHT");
    }
  case CEIL_2D:
  case EUC_2D:
  case ATT:
    {
      int scrap;
      double x, y;
      cities.reserve(prob_size);

      for (i=0; i<prob_size; i++) {
	in_file>>scrap>>x>>y;
	cities.push_back(pair<double,double>(x,y));  // load the cities coordinate
      }

      for (i=0; i<maxrow; i++) {
	d[i][i]=0;
	for (j=i+1; j<prob_size; j++) {
	  const double dx=cities[i].first  - cities[j].first;
	  const double dy=cities[i].second - cities[j].second;

	  switch (weight_type) {
	  case EUC_2D: {
	    const double euc_d = sqrt(dx*dx+dy*dy);
	    d[i][j] = (int)(euc_d+0.5);   // round to nearest integer
	    break;
	  }
	  case CEIL_2D: {
	    const double euc_d = sqrt(dx*dx+dy*dy);
	    d[i][j] = (int)(ceil(euc_d));  // CEIL_2D
	    break;
	  }
	  case ATT: {
	    const double euc_d = sqrt((dx*dx+dy*dy)/10.0);
	    const int tij = (int)(euc_d+0.5);
	    if (tij<euc_d) 
	      d[i][j]=tij+1;
	    else 
	      d[i][j] = tij;
	    break;
	  }
	  default:
	    break;
	  }
	  if (j<maxrow) d[j][i] = d[i][j];
	}
      }
      break; // case 1;
    }
  case EXPLICIT:
    { // EXPLICIT
      switch (weight_format) {
      case FORMAT_NONE: {
	cerr << "Invalid weight format" << endl;
	in_file.close();
	abort();
	//	throw("INVALID_WEIGHT");
      }
      case FULL_MATRIX: {
	// read explicit distance full matrix
	for (i=0; i<prob_size; i++) {
	  for (j=0; j<prob_size; j++) {
	    in_file >> d[i][j];
	  }
	}
	break;
      }
      case LOWER_DIAG_ROW: {
	int di;
	// read lower diag row matrix
	for (i=0; i<prob_size; i++) {
	  for (j=0; j<=i && j<prob_size; j++) {
	    in_file >> di;
	    d[i][j] =  di;
	    if (i!=j) d[j][i] = di;
	  }
	}
	break;
      }
      case UPPER_DIAG_ROW: {
	int di;
	// read upper diag row matrix
	for (i=0; i<prob_size; i++) {
	  for (j=i; j<prob_size; j++) {
	    in_file >> di;
	    d[i][j] =  di;
	    if (i!=j) d[j][i] = di;
	  }
	}
	break;
      }
      case UPPER_ROW: {
	int di;
	// read upper row matrix
	for (i=0; i<prob_size; i++) {
	  d[i][i]=0;
	  for (j=i+1; j<prob_size; j++) {
	    in_file >> di;
	    d[i][j] =  di;
	    d[j][i] = di;
	  }
	}
	break;
      }

      }
    }
  }
  in_file.close();
  cout << "File loaded" << endl;
  candidate.init(40);
}

// destructor. Free allocated memory for the distance matrix
tsp_prob::~tsp_prob()
{
  if (d!=0) {

    unsigned i;
    for (i=0; i<prob_size; i++)
      if (d[i])
	delete[] d[i];

    delete[] d;
  }
}


//Function to evaluate the solution
//according to the total distance of the solution
int tsp_prob::evaluation(const tour& sol) const
{
  int eval=0;
  for (unsigned i=1; i<size(); i++)
    eval+=dist(sol[i-1],sol[i]);

  eval+=dist(sol[sol.size()-1],sol[0]);

  return eval;
}


bool tsp_prob::is_valid(const tour& sol) const
{
  if (sol.size() != prob_size) 
    {
      cerr << "badsize:" <<sol.size()<< endl;
      return false;
    }

  vector<int> visite(sol.size(), 0);
  // toutes les villes doivent etre visite une seule fois.
  for (unsigned i=0; i!=prob_size; ++i)
    {
      if (sol[i]>=prob_size) {
	cerr << "out_of_range: " << sol[i] << endl;
	return false;
      }

      if (++visite[sol[i]] > 1) {
	cerr << "too_many" << endl;
	return false;
      }
      assert(sol.checkB());
    }

  return true;
}



void tsp_prob::canonical_sol(tour& sol) const
{
  vector<unsigned> tmp_sol;
  tmp_sol.reserve(prob_size+1);
  unsigned i=0;

  for (i=0; i<prob_size; i++)
    tmp_sol.push_back(i);

  sol = tmp_sol;
}


void tsp_prob::plot_sol(const tour& sol, ostream& x) const
{
  if (weight_type==EUC_2D || weight_type==CEIL_2D) {
    x << "# distance: " << evaluation(sol) << endl;
    for (unsigned i=0; i!=prob_size; ++i) {
      x << cities[sol[i]].first << " " << cities[sol[i]].second << endl;
    }
    x << cities[sol[0]].first << " " << cities[sol[0]].second << endl;

  }
}

