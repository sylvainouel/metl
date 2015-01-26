#include "tsp_prob.hh"

#include <set>
#include <map>
#include <deque>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>


#include "sprng_rand.hh"

#include "tsp_eax.hh"

using namespace std;

void gready_subtours_recombine(vector<deque<int> >& subtours);


struct edge {
  edge(int _v1=0, int _v2=0, int _parent=0) : v1(_v1), v2(_v2), parent(_parent) { }
    
  bool operator<(const edge& e) const {
    if (max(v1, v2) < max(e.v1, e.v2)) return true;
    if (max(v1, v2) > max(e.v1, e.v2)) return false;
      
    if (min(v1, v2) < min(e.v1, e.v2)) return true;
    if (min(v1, v2) > min(e.v1, e.v2)) return false;
      
    return (parent<e.parent);
  }
    
  bool operator==(const edge& e) const {
    if ((v1 == e.v1) && (v2==e.v2) && (parent==e.parent)) return true;
    if ((v1 == e.v2) && (v2==e.v1) && (parent==e.parent)) return true;
    return false;  
  }
    
  int v1;
  int v2;
  int parent;
};





void tsp_eax::operator()(const tsp_prob::soleval_t& Ap,
			 const tsp_prob::soleval_t& Bp,
			 tsp_prob::soleval_t& Cp) const
{
  const tsp_prob::sol_t& A = Ap.first;
  const tsp_prob::sol_t& B = Bp.first;
  tsp_prob::sol_t& C = Cp.first;

  // cout <<A.size()<<"   "<<B.size()<<endl;
  assert(A.size()==B.size());
  assert(tsp_prob::instance().is_valid(A));

  // associate (parent, vertex) to all other vertex connected by edges
  typedef pair<int,int> keytype;
  typedef multimap<keytype, int> Rtype;
  Rtype R;

  set<edge> Cedges;
    
  for (unsigned i=0; i<A.size()-1; ++i) {
    // this structure is used for construction of the AB-cycles
    // It is indexed so that it is easy to find all edges from parent
    // X and starting from vectex V
    R.insert(pair<keytype, int>(keytype(0,A.get_tour()[i]), A.get_tour()[i+1]));
    R.insert(pair<keytype, int>(keytype(0,A.get_tour()[i+1]), A.get_tour()[i]));
    R.insert(pair<keytype, int>(keytype(1,B.get_tour()[i]), B.get_tour()[i+1]));
    R.insert(pair<keytype, int>(keytype(1,B.get_tour()[i+1]), B.get_tour()[i]));

    Cedges.insert(edge(A.get_tour()[i],   A.get_tour()[i+1]));   // C starts with a copy of A
  }
  R.insert(pair<keytype, int>(keytype(0,A.get_tour()[A.size()-1]), A.get_tour()[0]));
  R.insert(pair<keytype, int>(keytype(1,B.get_tour()[B.size()-1]), B.get_tour()[0]));
  R.insert(pair<keytype, int>(keytype(0,A.get_tour()[0]), A.get_tour()[A.size()-1]));
  R.insert(pair<keytype, int>(keytype(1,B.get_tour()[0]), B.get_tour()[B.size()-1]));

  Cedges.insert(edge(A.get_tour()[A.size()-1],   A.get_tour()[0])); 


  int parent;
  set<edge> edge_in_subtour;   // both edge_in_subtour and
  vector<edge> edge_vect;      // edge_vect keep the same
  // thing. edge_in_subtour is indexed
  // while edge_vect is ordered

  // associate a city with each time it was visited in edge_vect
  multimap<int, int> visited_cities;
  int cycle=-1;
  int start_vertex;

  vector<vector<edge> > ABcycles;  // this is the list of usefull ABcycles

  while(R.size()>0) {
    visited_cities.clear();
    edge_in_subtour.clear();
    edge_vect.clear();

    Rtype::iterator start_edge = R.begin();
    advance(start_edge, metl::rng(R.size()-1));      // choisi une arrete aleatoirement dans R

    start_vertex = start_edge->first.second;     // ville de depart
    parent = start_edge->first.first;            // parent de cette arrete
    visited_cities.insert(pair<int,int>(start_vertex,0)); // 0 is edge_vect_size()

    vector<edge> available_edge;
    do {
      pair<Rtype::iterator, Rtype::iterator> range;
      range = R.equal_range(keytype(parent, start_vertex));       // trouve toutes les arretes de ce parent et partant de cette vertex

      available_edge.clear();

      for (Rtype::iterator it=range.first; it!=range.second; ++it) {
	const edge new_edge(start_vertex, it->second, parent);
	// remove used edges
	if (edge_in_subtour.find(new_edge)==edge_in_subtour.end()) {
	  available_edge.push_back(new_edge);
	}
      }

      assert(available_edge.size()!=0);

      int choice = metl::rng(available_edge.size());
      // this is my chosen edge
      const edge& chosen_edge=*(available_edge.begin() + choice);

      edge_in_subtour.insert(chosen_edge);
      edge_vect.push_back(chosen_edge);
      start_vertex = chosen_edge.v2;
    
      cycle=-1;
      if (visited_cities.count(chosen_edge.v2)>0) {
	pair<multimap<int,int>::iterator, multimap<int,int>::iterator> city_it = visited_cities.equal_range(chosen_edge.v2);
	for (multimap<int,int>::iterator it=city_it.first; it!=city_it.second; ++it)
	  if (edge_vect[it->second].parent==!parent) {
	    // found an AB-cycle
	    cycle = (*it).second;  // remember where the cycle begins
	    // in edge_vect so we can truncate
	    // everything else
	    break;
	  }
      }
      if (cycle==-1) {
	// marque la ville comme etant visitee a cette position dans le vecteur d'arrete
	visited_cities.insert(pair<int, int>(start_vertex,edge_vect.size()));
      }
      parent = !parent;   // change de parent
    } while (cycle==-1);   // jusqu'a ce qu'on trouve un cycle


    // troncate this cycle to keep only an AB-cycle

    int removed_count;
    // remove edges in the AB-cycle from R
    for (vector<edge>::iterator i=edge_vect.begin()+cycle; i!=edge_vect.end(); ++i)
      {
	removed_count=0;
	pair<Rtype::iterator, Rtype::iterator> range = R.equal_range(keytype((*i).parent, (*i).v1));
	for (Rtype::iterator it=range.first; it!=range.second; ++it) {
	  if ((*it).second == (*i).v2) {
	    R.erase(it);
	    removed_count++;
	    break;
	  }
	}
	range = R.equal_range(keytype((*i).parent, (*i).v2));
	for (Rtype::iterator it=range.first; it!=range.second; ++it) {
	  if ((*it).second == (*i).v1) {
	    R.erase(it);
	    removed_count++;
	    break;
	  }
	}
	//	cout << removed_count << endl;
	// 	cout << i->v1 << "  " << i->v2 << "  " << i->parent<< endl;
	assert(removed_count==2);
	//	 	cout << i->v1 << "  " << i->v2 << "  " << i->parent<< endl;
      }
    //      cout << "end"<< endl;
    if (edge_vect.end() - (edge_vect.begin()+cycle) >2) {
      // usefull AB-cycle
      ABcycles.push_back(vector<edge>(edge_vect.begin()+cycle, edge_vect.end()));
    }
  }
  
  //   cout << "number of usefull ABcycles: "<<  ABcycles.size() << endl;


  // Take a copy of A (Cedges) and process it with the edges in the E-set
  // Remove edges in the E-set from A and Add edges from B
  for (vector<vector<edge> >::iterator it=ABcycles.begin(); it!=ABcycles.end(); ++it) 
    {
      if (metl::rng()<0.5) 
	{     // consider ABcycle with proibability 0.5
	  for (vector<edge>::iterator j=it->begin(); j!=it->end(); ++j) 
	    {
	      if (j->parent==0) 
		{ // if parent of this edge is A
		  Cedges.erase(*j);
		} 
	      else 
		{
		  Cedges.insert(*j);
		}
	    }
	}
    }

  // find all subtours in C, output them in subtours, in lists of visited cities
  vector<int> visit(A.size(),-1);  // associate city with tour number
  vector<deque<int> > subtours;

  for (set<edge>::iterator itt=Cedges.begin(); itt!=Cedges.end(); ++itt) 
    {
      if (visit[itt->v1]==-1 && visit[itt->v2]==-1)
	{ 
	  // start a new subtour
	  visit[itt->v1] = visit[itt->v2] = subtours.size();
	  subtours.push_back(deque<int>());
	  subtours.back().push_back(itt->v1);
	  subtours.back().push_back(itt->v2);
	  continue;
	}
      if (visit[itt->v1]!=-1 && (visit[itt->v2]==visit[itt->v1])) 
	{
	  continue;
	}
      int v=max(visit[itt->v1], visit[itt->v2]);
    
      if (itt->v1 == subtours[v].front()) 
	{
	  // add at the begining
	  subtours[v].push_front(itt->v2);
	  if (visit[itt->v2]==-1) 
	    visit[itt->v2]=v;
	} 
      else
	if (itt->v2 == subtours[v].back()) 
	  {
	    // add at the end
	    subtours[v].push_back(itt->v1);
	    if (visit[itt->v1]==-1) 
	      visit[itt->v1]=v;
	  } 
	else
	  if (itt->v2 == subtours[v].front()) 
	    {
	      // add at the begining
	      subtours[v].push_front(itt->v1);
	      if (visit[itt->v1]==-1) 
		visit[itt->v1]=v;
	    }
	  else
	    if (itt->v1 == subtours[v].back()) 
	      {
		// add at the end
		subtours[v].push_back(itt->v2);
		if (visit[itt->v2]==-1) 
		  visit[itt->v2]=v;
	      }

      if (visit[itt->v1]==visit[itt->v2]) continue;  // do not merge with same tour

      // merge visit[i->v1] with visit[i->v2]
      int v1 = visit[itt->v1];
      int v2 = visit[itt->v2];
      if (subtours[v1].front() == subtours[v2].back())
	// add subtour[v1] after subtour[v2]
	for (deque<int>::iterator k=subtours[v1].begin()+1; k!=subtours[v1].end(); ++k) 
	  {
	    subtours[v2].push_back(*k);
	    visit[*k]=v2;
	  }
      else 
	if (subtours[v1].back() == subtours[v2].front())
	  // add subtour[v1] before subtour[v2]
	  for (deque<int>::reverse_iterator k=subtours[v1].rbegin()+1; k!=subtours[v1].rend(); ++k) 
	    {
	      subtours[v2].push_front(*k);
	      visit[*k]=v2;
	    } 
	else 
	  if (subtours[v1].front() == subtours[v2].front())
	    // add subtour[v1] before subtour[v2] reversing it
	    for (deque<int>::iterator k=subtours[v1].begin()+1; k!=subtours[v1].end(); ++k) 
	      {
		subtours[v2].push_front(*k);
		visit[*k]=v2;
	      }
	  else 
	    if (subtours[v1].back() == subtours[v2].back())
	      // add subtour[v1] after subtour[v2] reversing it
	      for (deque<int>::reverse_iterator k=subtours[v1].rbegin()+1; k!=subtours[v1].rend(); ++k) 
		{
		  subtours[v2].push_back(*k);
		  visit[*k]=v2;
		}
      subtours[v1].clear();  //erase v1
    }

  // remove empty subtours  
  for (vector<deque<int> >::iterator ii=subtours.begin(); ii!=subtours.end();) 
    {
      if (ii->empty()) 
	{
	  subtours.erase(ii);
	  continue;
	}
      ++ii;
    }

  //  cout << "number of subtours: " << subtours.size() << endl; 
  if (subtours.size()>1) 
    {
      gready_subtours_recombine(subtours);
    }
  //  C.clear();  // C is the output vector
  vector<unsigned> tmp_sol;
  tmp_sol.reserve(A.size());

  for(deque<int>::iterator id=subtours.front().begin(); id!=subtours.front().end(); ++id)
    {
      tmp_sol.push_back(*id);
    }


  C = tour(tmp_sol);
  Cp.second = tsp_prob::instance().evaluation(Cp.first);
}



void tsp_eax::gready_subtours_recombine(vector<deque<int> >& subtours) const {

  int mincost;
  while (subtours.size()>1) {
    mincost=std::numeric_limits<int>::max();
    vector<deque<int> >::iterator other=subtours.begin();
    int other_pos = 0;
    int shortest_pos = 0;


    //    vector<deque<int> >::iterator shortest=subtours.begin()+metl::rng(subtours.size());
    //     // find shortest subtours (in terms of number of edges)
    vector<deque<int> >::iterator shortest=subtours.begin();

    for (vector<deque<int> >::iterator it=subtours.begin()+1; it!=subtours.end(); ++it) {
      if (shortest->size() > it->size()) shortest=it;
    }
    

    unsigned size1 = shortest->size();
    int type=0;

    const tsp_prob& problem = tsp_prob::instance();

    for (unsigned i=0; i<size1; ++i) { // consider each edge of shortest

      //       vector<deque<int> >::iterator it=subtours.begin()+metl::rng(subtours.size());
      //       while(it==shortest) it=subtours.begin()+metl::rng(subtours.size());

      // for each other tour
      for (vector<deque<int> >::iterator it=subtours.begin(); it!=subtours.end(); ++it) {
	if (it==shortest) continue;

	// consider each edge of other subtour
	for (unsigned j=0; j<it->size(); ++j) {
	  const int cost_base = 
	    -problem.dist(it->front(), it->back()) 
	    -problem.dist(shortest->front(), shortest->back());

	  const int cost =
	    cost_base+
	    problem.dist(it->front(), shortest->back()) +
	    problem.dist(it->back(), shortest->front());
	  const int cost2 = 
	    cost_base+
	    problem.dist(it->front(), shortest->front()) +
	    problem.dist(it->back(), shortest->back());
	    
	  //	  cout << "cost: " << cost << " "<< it->front() << " " << it->back() <<"  " << shortest->front() << " " << shortest->back() <<endl;
	  if (cost < mincost) {
	    other = it;
	    other_pos = j;
	    shortest_pos = i;
	    mincost = cost;
	    type=0;
	  }
	  if (cost2 < mincost) {
	    other = it;
	    other_pos = j;
	    shortest_pos = i;
	    mincost = cost2;
	    type=1;
	  }


	  it->push_back(it->front());
	  it->pop_front();
	}
      }

      shortest->push_back(shortest->front());
      shortest->pop_front();
    }

    assert(shortest!=other);
    //    cout << "merge: " << shortest - subtours.begin() << " with: " << other - subtours.begin() << endl;

    // merge the two subtours
    for (int i=0; i<shortest_pos; ++i)
      {
	shortest->push_back(shortest->front());
	shortest->pop_front();
      }

    for (int j=0; j<other_pos; ++j)
      {
	other->push_back(other->front());
	other->pop_front();
      }
    if (type==1) {
      reverse(other->begin(), other->end());
    }

    // copy shortest after other
    for (deque<int>::iterator idd=shortest->begin(); idd!=shortest->end(); ++idd) {
      other->push_back(*idd);
    }
    subtours.erase(shortest);
    //    cout << "subtours size: " << subtours.size() << endl;
  }

  //  cout << "final tour:";
  //   for (deque<int>::iterator i=subtours[0].begin(); i!=subtours[0].end(); ++i) {
  //     cout << *i << " ";
  //   }
  //   cout << endl;

}


